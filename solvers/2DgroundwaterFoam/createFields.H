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

Info << "Computing (or reading if present) hwater" << endl;
volScalarField hwater
(
    IOobject
    (
        "hwater",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    potential
);
hwater == potential - z0;

if (min(hwater).value() < 0)
{
    Info << "*** WARNING *** : hwater negative (initialized to 0.01) " << endl;

    forAll(mesh.C(),celli)
    {
        hwater[celli] = max(hwater[celli],0.01);
    }
    forAll(hwater.boundaryField(),patchi)
    {
        forAll(hwater.boundaryFieldRef()[patchi],facei)
        {
            hwater.boundaryFieldRef()[patchi][facei] =  max(0.01,hwater.boundaryFieldRef()[patchi][facei]);
        }
    }
}
hwater.write();

Info << "min(hwater) = " << min(hwater).value() << " ; max(hwater) = " << max(hwater).value() << nl << endl;

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

const dimensionedScalar nu(transportProperties.lookup("nu"));

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
Info<< "Reading field K\n" << endl;
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
volScalarField M(K*g/nu);
surfaceScalarField Kf = fvc::interpolate(K,"K");
surfaceScalarField Mf("Mf",Kf*g/nu);

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

// field necessary for darcyGradPressure boundary condition
surfaceScalarField phiG("phiG",0*phi);
surfaceScalarField phiPc("phiPc",0*phi);

//-
Info << "Reading field infiltration (if present)\n" << endl;
dimensionedScalar infiltrationScalar(transportProperties.lookupOrDefault("infiltration",dimensionedScalar("",dimLength/dimTime,0.)));
volScalarField infiltration
(
    IOobject
    (
        "infiltration",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    infiltrationScalar
);
