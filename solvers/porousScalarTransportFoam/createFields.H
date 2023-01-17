////////////////////////////////////////////////////////////////
////////////////////// GENERAL PROPERTIES //////////////////////
////////////////////////////////////////////////////////////////

Info << "Reading transportProperties" << endl;
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

//- list that receives event files for active event-based boundary conditions
List<patchEventFile*> patchEventList;
eventFlux::setEventFileRegistry(&patchEventList, "C");

List<sourceEventFile*> sourceEventList;

wordList speciesNames(transportProperties.lookupOrDefault("species", wordList(1, "C")));

Info<< "Reading composition\n" << endl;

multiscalarMixture composition
(
    transportProperties,
    speciesNames,
    mesh,
    word::null,
    eps,
    &sourceEventList,
    "C"
);

///////////////////////////////////////////////////////////////////
///////////////////// VELOCITY - FLUXES ///////////////////////////
///////////////////////////////////////////////////////////////////

Info<< "Creating fluidPhase" << endl;
word phaseName(transportProperties.optionalSubDict("porousTransport").lookupOrDefault<word>("phaseName",""));
autoPtr<fluidPhase> phase = fluidPhase::New(mesh, transportProperties, phaseName);
volVectorField& U = phase->U();
surfaceScalarField& phi = phase->phi();

///////////////////////////////////////////////////////////////////
////////////////////// POTENTIAL-SATURATION ///////////////////////
///////////////////////////////////////////////////////////////////

Info << nl << "Reading water content theta and/or Saturation field S" << phaseName << "..." << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("theta",dimless,1)
);

volScalarField Saturation
(
    IOobject
    (
        "S"+phaseName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("S"+phaseName,dimless,1)
);

word Uname="U"+phaseName;
if (theta.headerOk())
{
    Saturation = theta/eps;
    Uname = "U";
    Info << "===> field S" << phaseName << " computed from theta/eps, min = " << gMin(Saturation.internalField()) << " ; max = " << gMax(Saturation.internalField()) << endl;
}
else
{
    if (Saturation.headerOk())
    {
        Info << "===> field S" << phaseName << " read, min = " << gMin(Saturation.internalField()) << " ; max = " << gMax(Saturation.internalField()) <<  endl;
    }
    else
    {
        Info << "===> Saturation file not found (saturated flow)" << endl;
    }
}

////////////////////////////////////////////////////
//////////////////// OUTPUT CSV ////////////////////
////////////////////////////////////////////////////

bool CSVoutput=runTime.controlDict().lookupOrDefault<bool>("CSVoutput",true);
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
