///////////////////////////////////////////////////////////////
///////////////////// POTENTIAL READING ///////////////////////
///////////////////////////////////////////////////////////////

Info << "Reading Pressure head h" << endl;
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

Info << "Creating deltah field (for Newton loop)" << endl;
volScalarField deltah
(
    IOobject
    (
        "deltah",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    h
);
deltah == dimensionedScalar("",dimLength,0);

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
eventInfiltration::setEventFileRegistry(&patchEventList, h.name());

/////////////////////////////////////////////////////////////////////////////
//////////////// RELATIVE PERMABILITY - CAPILLARY MODEL /////////////////////
/////////////////////////////////////////////////////////////////////////////

//- theta temporary initialisation
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
    dimless,
    calculatedFvPatchScalarField::typeName
);
theta.dimensions().reset(dimless);
theta = 0;

// Relative permeability / Capillary pressure models
autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New("krModel",transportProperties,theta);
autoPtr<capillarityModel> pcModel = capillarityModel::New("pcModel",transportProperties,theta);

Info << nl << "Computing saturation field theta" << endl;
theta = pcModel->correctAndSb(h);
theta.write();

////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL ////////////////////////////
////////////////////////////////////////////////////////////////////

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
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Intrinsic permeability
Info  << nl << "Reading permeability field K" << endl;
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

// permeability interpolation
surfaceScalarField Kf(fvc::interpolate(K,"K"));

// specific storage
dimensionedScalar Ss(transportProperties.lookupOrDefault<dimensionedScalar>("Ss",dimensionedScalar("Ss",dimless/dimLength,0.)));
Info << nl << "Reading speficic storage : Ss = " << Ss.value() << endl;

///////////////////////////////////////////////////////////////////
////////////////////////// FORCING TERMS //////////////////////////
///////////////////////////////////////////////////////////////////

//- For constant or event source injection/extraction
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

////////////////////////////////////////////////////
//////////////////// OUTPUT CSV ////////////////////
////////////////////////////////////////////////////

bool CSVoutput=runTime.controlDict().lookupOrDefault<bool>("CSVoutput",true);
OFstream waterMassBalanceCSV("waterMassBalance.csv");
if (CSVoutput)
{
    waterMassBalanceCSV << "#Time ";
    forAll(mesh.boundaryMesh(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            waterMassBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
        }
    }
    waterMassBalanceCSV << endl;
}
