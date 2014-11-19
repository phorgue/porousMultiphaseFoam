//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

bool CFL =  runTime.controlDict().lookupOrDefault<bool>("CFL", true);

Info<< "Reading field p\n" << endl;
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



Info<< "Reading field Sb \n" << endl;
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

Info<< "Reading transportProperties\n" << endl;

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

/////////////////////////////////////////////////////////////////////////////
///////////////////////// SATURATION INITIALIZATION /////////////////////////
/////////////////////////////////////////////////////////////////////////////

//Minimal and Maximal saturation (physical condition)
dimensionedScalar Smin(transportProperties.lookup("Smin"));
dimensionedScalar Smax(transportProperties.lookup("Smax"));
	
// Test  the initial Saturation if inferior to Smin or superior to Smax.
if (min(Sb).value() < Smin.value()) 
{
    FatalErrorIn
        (
            "in createFields.H"
        )
        << "Saturation equal or inferior to minimal saturation " 
            << exit(FatalError);
}
if ( max(Sb).value() > Smax.value()) 
{  
    FatalErrorIn
        (
            "in createFields.H"
        )
        << "Saturation equal or superior to maximal saturation " 
            << exit(FatalError);
}

forAll(Sb,celli)
{
    if(Sb[celli] == Smin.value()) Sb[celli] += SMALL;
    if(Sb[celli] == Smax.value()) Sb[celli] -= SMALL;
}

// Relative permeability model 
autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New("krModel",transportProperties,Sb);

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Reading Porous Media Properties \n" << endl;
// Porosity	
dimensionedScalar eps(transportProperties.lookup("eps"));
   
// Intrinsic permeability       
Info<< "Reading field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// permeability interpolation (harmonic by default)
surfaceScalarField Kf = fvc::interpolate(K,"K");

/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Reading field U\n" << endl;
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

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// CAPILLARY MODEL /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));
autoPtr<capillarityModel> myCapModel = capillarityModel::New("pc",transportProperties,Sb);

volScalarField cellVolume
(
    IOobject
    (
        "cellVolume",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("unite",dimVolume,1.0)//pour la division dans les termes sources
);
cellVolume.internalField() = mesh.V();

