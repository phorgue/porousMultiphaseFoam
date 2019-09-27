//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info << "Reading transportProperties" << endl;
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

//- dictionary for handling event time from boundary conditions
IOdictionary patchEventFiles
(
    IOobject
    (
        "patchEventFiles",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar rho(transportProperties.lookup("rho"));
const dimensionedScalar mu(transportProperties.lookup("mu"));

/////////////////////////////////////////////////////////////////
////////////////////// HEIGHT - POTENTIAL ///////////////////////
/////////////////////////////////////////////////////////////////

Info << nl << "Reading field potential" << endl;
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

Info << nl << "Reading elevation field z0" << endl;
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
Info << "min(z0) = " << min(z0).value() << " ; max(z0) = " << max(z0).value() << endl;

Info << nl << "Computing hwater" << endl;
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
<< " ; hwaterMin = " << hwaterMin.value() << endl;

hwater.write();

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Porosity
Info<< "Reading porosity field eps (if present)" << endl;
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
Info << nl << "Reading permeability field K" << endl;
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

Info<< "Reading field U" << nl << endl;
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
surfaceScalarField phiWater = phi * fvc::interpolate(hwater);

// field necessary for darcyGradPressure boundary condition
surfaceScalarField phiG("phiG",0*phi);
surfaceScalarField phiPc("phiPc",0*phi);

/////////////////////////////////////////////////////////////////////////
////////////////////////// FORCING TERMS ////////////////////////////////
/////////////////////////////////////////////////////////////////////////

Info << "Reading infiltration field (if present)" << endl;
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

//- Tracer injection source term
Info << "Reading tracer source field (if present)" << endl;
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
    dimensionedScalar("",dimMass*dimLength/(dimVol*dimTime),0)
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

///////////////////////////////////////////////////////////////////////
////////////////////////// TRACER PART ////////////////////////////////
///////////////////////////////////////////////////////////////////////

Info << nl << "Reading concentration field C" << endl;
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

Info << nl << "Reading porous transport properties" << endl;
dictionary porousTransport(transportProperties.subDict("porousTransport"));
word phaseName(porousTransport.lookupOrDefault<word>("phaseName",""));
dimensionedScalar Kd(porousTransport.lookupOrDefault("Kd",dimensionedScalar("Kd_default",dimensionSet(-1,3,0,0,0),0.)));
dimensionedScalar rs(porousTransport.lookupOrDefault("rs",dimensionedScalar("rs_default",dimensionSet(1,-3,0,0,0),0.)));
dimensionedScalar epsTotal(porousTransport.lookupOrDefault("epsTotal",dimensionedScalar("epsTotal_default",dimless,0.)));
dimensionedScalar lambda(porousTransport.lookupOrDefault("lambda",dimensionedScalar("lambda_default",dimensionSet(0,0,-1,0,
0),0.)));

Info << nl << "porousTransport parameters" << nl << "{" << endl;
Info << "    " << "phaseName " << " : " << phaseName << endl;
Info << "    " << Kd.name() << " : " << Kd.value() << endl;
Info << "    " << rs.name() << " : " << rs.value() << endl;
Info << "    " << epsTotal.name() << " : " << epsTotal.value() << endl;
Info << "    " << lambda.name() << " : " << lambda.value() << endl;
Info << "}" << endl;

//- creation of dispersion model
autoPtr<dispersionModel> DeffModel= dispersionModel::New("DeffModel",transportProperties, mesh);
volTensorField Deff =  eps*hwater*DeffModel->Deff();

Info << nl << "Computing retard coefficient R = ";
volScalarField R =  1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps;
Info << R << endl;

////////////////////////////////////////////////////
//////////////////// OUTPUT CSV ////////////////////
////////////////////////////////////////////////////

bool CSVoutput=runTime.controlDict().lookupOrDefault<bool>("CSVoutput",true);
OFstream waterMassBalanceCSV("waterMassBalance.csv");
if (CSVoutput)
{
    waterMassBalanceCSV << "#Time flux(Infiltration) flux(fixedPoints)";
    forAll(mesh.boundaryMesh(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            waterMassBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
        }
    }
    waterMassBalanceCSV << endl;
}

OFstream CmassBalanceCSV("CmassBalance.csv");
if (CSVoutput)
{
    CmassBalanceCSV << "#Time TotalMass(kg)";
    forAll(mesh.boundaryMesh(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            CmassBalanceCSV << " flux(" << phiWater.boundaryField()[patchi].patch().name() << ")";
        }
    }
    CmassBalanceCSV << " flux(fixedPoint)" << endl;
}
