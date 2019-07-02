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

const dimensionedScalar rho(transportProperties.lookup("rho"));
const dimensionedScalar mu(transportProperties.lookup("mu"));

//////////////////////////////////////////////////////////////////
////////////////////// Height - Potential ///////////////////////
//////////////////////////////////////////////////////////////////

Info<< "Reading field potential\n" << endl;
volScalarField potential
(
    IOobject
    (
        "potential",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field z0" << endl;
volScalarField z0
(
    IOobject
    (
        "z0",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info << "min(z0) = " << min(z0).value() << " ; max(z0) = " << max(z0).value() << nl << endl;

Info << "Computing hwater" << endl;
volScalarField hwater
(
    IOobject
    (
        "hwater",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    potential
);
hwater == potential - z0;

dimensionedScalar hwaterMin(transportProperties.lookupOrDefault("hwaterMin",dimensionedScalar("hwaterMin",dimLength,0.01)));
scalar cumulativeWaterAdded = 0;

#include "correctInitialPotential.H"

Info << "min(hwater) = " << min(hwater).value() << " ; max(hwater) = " << max(hwater).value()
<< " ; hwaterMin = " << hwaterMin.value() << nl << endl;

hwater.write();

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
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    epsScalar
);

if (max(eps).value() <= 0) FatalErrorIn("createFields.H") << " porosity eps equal to 0" << exit(FatalError);

// Intrinsic permeability
Info<< "Reading field K" << endl;
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
Info << "min(K) = " << min(K).value() << " ; max(K) = " << max(K).value() << nl << endl;

// permeability/mobility field
dimensionedScalar g("g",dimLength/(dimTime*dimTime),9.81);
volScalarField M(K*rho*g/mu);
surfaceScalarField Kf = fvc::interpolate(K,"K");
surfaceScalarField Mf("Mf",Kf*g*rho/mu);

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
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"
phi == (-Mf * fvc::snGrad(potential)) * mesh.magSf();
forAll(mesh.boundary(),patchi)
{
    if (isA< fixedValueFvPatchField<vector> >(U.boundaryField()[patchi]))
    {
        phi.boundaryFieldRef()[patchi] = U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    }
}
surfaceScalarField transmissivity("transmissivity",Mf*fvc::interpolate(hwater));

// field necessary for darcyGradPressure boundary condition
surfaceScalarField phiG("phiG",0*phi);
surfaceScalarField phiPc("phiPc",0*phi);

//-
Info << "Reading field infiltration (if present)" << endl;
dimensionedScalar infiltrationScalar(transportProperties.lookupOrDefault("infiltration",dimensionedScalar("",dimLength/dimTime,0.)));
volScalarField infiltration
(
    IOobject
    (
        "infiltration",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    infiltrationScalar
);

volScalarField seepageTerm
(
    IOobject
    (
        "seepageTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("seepage_value",dimLength/dimTime,0.)
);

scalar zScale(mesh.bounds().max().z()-mesh.bounds().min().z());

//- Looking for 2nd order time-scheme for source term
bool isBackwardScheme = false;
if (word(mesh.ddtScheme("source")) == "backward")
{
    isBackwardScheme = true;
}

////////////////////////////////////////////////////
//////////////////// OUTPUT CSV ////////////////////
////////////////////////////////////////////////////

bool CSVactive=runTime.controlDict().lookupOrDefault<bool>("CSVactive",true);
word CSVoutput=runTime.controlDict().lookupOrDefault<word>("CSVoutput","waterMassBalance.csv");
OFstream massBalanceCSV(CSVoutput);
if (CSVactive)
{
    massBalanceCSV << "#Time flux(Infiltration) flux(fixedPoints)";
    forAll(mesh.boundaryMesh(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            massBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
        }
    }
    massBalanceCSV << endl;
}
