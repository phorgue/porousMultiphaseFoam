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

///////////////////////////////////////////////////////////////////////
////////////////////////// SCALAR TRANSPORT  //////////////////////////
///////////////////////////////////////////////////////////////////////

Info<< "Reading porosity eps\n" << endl;
dimensionedScalar eps
(
    transportProperties.lookupOrDefault("eps",dimensionedScalar("eps_default",dimless,1))
);

Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading porous transport properties" << endl;
dictionary porousTransport(transportProperties.subDict("porousTransport"));
word phaseName(porousTransport.lookupOrDefault<word>("phaseName","notDefined"));
dimensionedScalar Kd(porousTransport.lookupOrDefault("Kd",dimensionedScalar("Kd_default",dimensionSet(-1,3,0,0,0),0.)));
dimensionedScalar rs(porousTransport.lookupOrDefault("rs",dimensionedScalar("rs_default",dimensionSet(1,-3,0,0,0),0.)));
dimensionedScalar epsTotal(porousTransport.lookupOrDefault("epsTotal",dimensionedScalar("epsTotal_default",dimless,0.)));
dimensionedScalar lambda(porousTransport.lookupOrDefault("lambda",dimensionedScalar("lambda_default",dimensionSet(0,0,-1,0,0),0.)));

Info << nl << "porousTransport parameters" << nl << "{" << endl;
Info << "    " << "phaseName " << " : " << phaseName << endl;
Info << "    " << Kd.name() << " : " << Kd.value() << endl;
Info << "    " << rs.name() << " : " << rs.value() << endl;
Info << "    " << epsTotal.name() << " : " << epsTotal.value() << endl;
Info << "    " << lambda.name() << " : " << lambda.value() << endl;
Info << "} \n" << endl;

//- creation of dispersion model
autoPtr<dispersionModel> DeffModel= dispersionModel::New("DeffModel",transportProperties,Utheta);

Info << "Computing retard coefficient R" << endl;
scalar R =  1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps.value();

#include "createWellbores.H"
dimensionedScalar Cin(wellboreProperties.lookup("Cin"));
