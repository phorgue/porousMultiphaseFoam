///////////////////////////////////////////////////////////////
////////////////////// TIMESTEP CONTROL ///////////////////////
///////////////////////////////////////////////////////////////

word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Coats");
if (adjustTimeStep) adjustTimeStep=true; // to remove warnings at compilation
int CFLint = -1;
if (CFL == "Coats") CFLint = 0;
else if (CFL == "Courant") CFLint = 1; 
else
{
    FatalErrorIn
        (
            "in createFields.H"
        )
        << "CFL condition unrecongnized : Coats and Courant available" 
            << exit(FatalError);
}

//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

Info << nl << "Reading field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << nl << "Reading field Sb" << endl;
volScalarField Sb
(
    IOobject
    (
        "Sb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax",0.));

//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info << nl << "Reading transportProperties" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//- list that receives event files of event-based boundary conditions
List<patchEventFile*> patchEventList;
eventInfiltration::setEventFileRegistry(&patchEventList, Sb.name());

dimensionedScalar Sbmin(transportProperties.lookupOrDefault("Sbmin",dimensionedScalar("Sbmin",dimless,0)));

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////

autoPtr<incompressiblePhase> phasea = incompressiblePhase::New(mesh,transportProperties,"a");
volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& mua = phasea->mu();

autoPtr<incompressiblePhase> phaseb = incompressiblePhase::New(mesh,transportProperties,"b");
volVectorField& Ub = phaseb->U();
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& mub = phaseb->mu();	

dimensionedScalar Mmu(mub/mua);//viscosity ratio
 
// Relative permeability / Capillary pressure models 
autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New(mesh, transportProperties, "Sb");
scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));
autoPtr<capillarityModel> pcModel = capillarityModel::New(mesh, transportProperties, "Sb");

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Porosity
Info << nl << "Reading porosity field eps (if present)" << endl;
dimensionedScalar epsScalar(transportProperties.lookupOrDefault("eps",dimensionedScalar("",dimless,0.)));
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    epsScalar
);

// Intrinsic permeability       
Info << nl << "Reading permeability field K" << endl;
volTensorField K
(
    IOobject
    (
        "K",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// interpolated permeability
surfaceTensorField Kf(fvc::interpolate(K,"K"));

/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info << nl << "Reading field U" << nl << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Ua + Ub
);

#include "createPhi.H"
surfaceScalarField phiP = phi;

///////////////////////////////////////////////////////////////////
////////////////////////// FORCING TERMS //////////////////////////
///////////////////////////////////////////////////////////////////

volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimless/dimTime,0)
);
