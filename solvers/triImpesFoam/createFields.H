//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Coats");
if (adjustTimeStep) adjustTimeStep=true; // to remove warnings at compilation  

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

Info<< "Reading field Sw \n" << endl;
volScalarField Sw
(
    IOobject
    (
        "Sw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field So \n" << endl;
volScalarField So
(
    IOobject
    (
        "So",
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

dimensionedScalar Swr(transportProperties.lookupOrDefault("Swr",dimensionedScalar("Swr",dimless,0)));

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////

autoPtr<incompressiblePhase> phasea = incompressiblePhase::New(mesh,transportProperties,"a");
volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& mua = phasea->mu();

autoPtr<incompressiblePhase> phaseo = incompressiblePhase::New(mesh,transportProperties,"o");
volVectorField& Uo = phaseo->U();
surfaceScalarField& phio = phaseo->phi();
const dimensionedScalar& rhoo = phaseo->rho();
const dimensionedScalar& muo = phaseo->mu();    

autoPtr<incompressiblePhase> phasew = incompressiblePhase::New(mesh,transportProperties,"w");
volVectorField& Uw = phasew->U();
surfaceScalarField& phiw = phasew->phi();
const dimensionedScalar& rhow = phasew->rho();
const dimensionedScalar& muw = phasew->mu();    

dimensionedScalar Mmu(muw/mua);//viscosity ratio

// Relative permeability model 
autoPtr<triRelativePermeabilityModel> krModel = triRelativePermeabilityModel::New("krModel",transportProperties,Sw, So);

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Reading Porous Media Properties \n" << endl;

// Porosity	
dimensionedScalar epsScalar(transportProperties.lookupOrDefault("eps",dimensionedScalar("",dimless,0.)));
   
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    epsScalar
);

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
    Ua + Uw + Uo
);

#include "createPhi.H"
surfaceScalarField phiP = phi;

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// CAPILLARY MODEL /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));

autoPtr<triCapillarityModel> pcModel = triCapillarityModel::New("pcModel",transportProperties,Sw,So);
