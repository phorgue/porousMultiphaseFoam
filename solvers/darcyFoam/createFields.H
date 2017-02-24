Info << "Reading field p\n" << endl;
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

autoPtr<incompressiblePhase> phase = incompressiblePhase::New(mesh,transportProperties,"a");
volVectorField& Ua = phase->U();
const dimensionedScalar& rho = phase->rho();
const dimensionedScalar& mu = phase->mu();

volVectorField U("U",Ua);

#include "createPhi.H"

Info << "Reading field K\n" << endl;
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

surfaceScalarField Kf(fvc::interpolate(K,"K"));
surfaceScalarField Mf ("Mf",Kf/mu);
surfaceScalarField phiG("phiG",(rho * Mf * g) & mesh.Sf());
surfaceScalarField phiPc("phiPc", 0 * phiG); //- required for porousBoundaryConditions
