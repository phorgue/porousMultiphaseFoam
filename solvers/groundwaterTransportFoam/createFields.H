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

Info<< "Reading transportProperties" << endl;
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
eventFlux::setEventFileRegistry(&patchEventList, "C");

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
//////////////// RELATIVE PERMABILITY - CAPILLARY MODEL /////////////////////
/////////////////////////////////////////////////////////////////////////////

// Relative permeability / Capillary pressure models
autoPtr<twophasePorousMediumModel> porousModel = twophasePorousMediumModel::New("theta", mesh, transportProperties, phasetheta);
autoPtr<relativePermeabilityModel>& krModel = porousModel->krModel();
autoPtr<capillarityModel>& pcModel = porousModel->pcModel();

// Water content computation
Info<< "Computing field theta \n" << endl;
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
theta = pcModel->correctAndSb(h);
theta.write();

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Intrinsic permeability
porousModel->check_K();
const volScalarField& K = porousModel->K();

// permeability interpolation (harmonic by default)
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

///////////////////////////////////////////////////////////////////////
////////////////////////// SCALAR TRANSPORT  //////////////////////////
///////////////////////////////////////////////////////////////////////

Info<< "Reading composition" << endl;
autoPtr<porousMediumTransportModel> pmTransportModel = porousMediumTransportModel::New("theta", porousModel);
multiscalarMixture& composition = pmTransportModel->composition();
List<sourceEventFile*>& tracerSourceEventList = pmTransportModel->sourceEventList();

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

PtrList<OFstream> CmassBalanceCSVs;
if (CSVoutput)
{
    CmassBalanceCSVs.resize(composition.Y().size());

    forAll(composition.Y(), speciesi)
    {
        CmassBalanceCSVs.set(speciesi, new OFstream(composition.species()[speciesi] + "massBalance.csv"));

        auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];

        CmassBalanceCSV << "#Time TotalMass(kg)";
        forAll(mesh.boundaryMesh(),patchi)
        {
            if (mesh.boundaryMesh()[patchi].type() == "patch")
            {
                CmassBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
            }
        }
        CmassBalanceCSV << endl;
    }
}

bool writeResiduals=false;
