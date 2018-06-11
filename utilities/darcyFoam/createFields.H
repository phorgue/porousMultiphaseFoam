word phaseName("a");
if (args.optionFound("phase")) phaseName = args.option("phase");

Info << nl << "Reading field p\n" << endl;
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

Info << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<incompressiblePhase> phase = incompressiblePhase::New(mesh,transportProperties,phaseName);
volVectorField& UphaseName = phase->U();
const dimensionedScalar& rho = phase->rho();
const dimensionedScalar& mu = phase->mu();
volVectorField U("U",UphaseName);

#include "createPhi.H"
surfaceScalarField phiPhaseName
(
    IOobject
    (
        "phi"+phaseName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
);

Info << "Reading field K\n" << endl;
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

surfaceScalarField Kf(fvc::interpolate(K,"K"));
surfaceScalarField Mf ("Mf",Kf/mu);
surfaceScalarField phiG("phiG",(rho * Mf * g) & mesh.Sf());
surfaceScalarField phiPc("phiPc", 0 * phiG); //- required for porousBoundaryConditions
