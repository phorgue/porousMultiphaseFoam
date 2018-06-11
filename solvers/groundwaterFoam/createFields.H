//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

Info << "Reading Pressure head h \n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

autoPtr<incompressiblePhase> phasetheta = incompressiblePhase::New(mesh,transportProperties,"theta");
volVectorField& Utheta = phasetheta->U();
const dimensionedScalar& rhotheta = phasetheta->rho();
const dimensionedScalar& mutheta = phasetheta->mu();    
phasetheta->phi().writeOpt()=IOobject::NO_WRITE;


surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Utheta) & mesh.Sf()
);

/////////////////////////////////////////////////////////////////////////////
///////////////////////// SATURATION INITIALIZATION /////////////////////////
/////////////////////////////////////////////////////////////////////////////

//Minimal and Maximal saturation (physical condition)
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimless
);
theta.dimensions().reset(dimless);
theta == 0; // temporary initialisation

/////////////////////////////////////////////////////////////////////////////
//////////////// RELATIVE PERMABILITY - CAPILLARY MODEL /////////////////////
/////////////////////////////////////////////////////////////////////////////

// Relative permeability model 
autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New("krModel",transportProperties,theta);
autoPtr<capillarityModel> pcModel = capillarityModel::New("pcModel",transportProperties,theta);

Info<< "Computing field theta \n" << endl;
theta == pcModel->correctAndSb(h);
theta.write();

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////
   
// Intrinsic permeability       
Info<< "Reading permeability field K\n" << endl;
volScalarField K
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

// permeability interpolation (harmonic by default)
surfaceScalarField Kf = fvc::interpolate(K,"K");
