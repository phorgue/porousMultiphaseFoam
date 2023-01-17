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

//- list that receives pointers to the event files of event-based boundary conditions
List<patchEventFile*> patchEventList;
eventFlux::setEventFileRegistry(&patchEventList, "C");

Info<< "Reading composition\n" << endl;
wordList speciesNames(transportProperties.lookupOrDefault("species", wordList(1, "C")));

List<sourceEventFile*> sourceEventList;
multiscalarMixture composition
(
    transportProperties,
    speciesNames,
    mesh,
    word::null,
    eps,
    &sourceEventList,
    "C",
    dimLength
);

/////////////////////////////////////////////////////////////////
////////////////////// HEIGHT - POTENTIAL ///////////////////////
/////////////////////////////////////////////////////////////////

Info << nl << "Reading hwater field (if present)" << endl;
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
    mesh,
    dimensionedScalar("hwater_default",dimLength,0.)
);

if (hwater.headerOk())
{
    Info << "min(hwater) = " << gMin(hwater.internalField()) << " ; max(hwater) = " << gMax(hwater.internalField()) << endl;
}
else
{
    Info << "file hwater not found " << endl
        << nl << "Reading potential field..." << endl;
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

    Info<< nl << "Reading field z0" << endl;
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
    Info << "min(z0) = " << gMin(z0.internalField()) << " ; max(z0) = " << gMax(z0.internalField()) << endl;

    hwater = potential - z0;

    Info << nl << "Computed hwater : min(hwater) = " << gMin(hwater.internalField()) << " ; max(hwater) = " << gMax(hwater.internalField()) << endl;
    hwater.write();
}

if (gMin(hwater.internalField()) < 0)
{
    FatalErrorIn("createFields") << " hwater has negative values" << abort(FatalError);
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Creating fluidPhase" << endl;
word phaseName(transportProperties.optionalSubDict("porousTransport").lookupOrDefault<word>("phaseName",""));
autoPtr<fluidPhase> phase = fluidPhase::New(mesh, transportProperties, phaseName);
volVectorField& U = phase->U();
surfaceScalarField& phi = phase->phi();
surfaceScalarField phihwater("phihwater", phi * fvc::interpolate(hwater));

///////////////////////////////////////////////////////////////////
////////////////////////// FORCING TERMS //////////////////////////
///////////////////////////////////////////////////////////////////

scalar zScale(mesh.bounds().max().z()-mesh.bounds().min().z());

Info << nl << "Read seepage flow rate if present...";
volScalarField seepageTerm
(
    IOobject
    (
        "seepageTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("seepage_value",dimLength/dimTime,0.)
);

if (seepageTerm.headerOk())
{
    Info << "OK" << endl;
    //- ensuring that there is no seepage inflow
    forAll(seepageTerm,celli)
    {
        seepageTerm[celli] = max(0,seepageTerm[celli]);
    }
}
else
{
    Info << " no seepageTerm file" << endl;
}

/////////////////////////////////////////////////////////
//////////////////// OUTLET CSV FILE ////////////////////
/////////////////////////////////////////////////////////

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
                CmassBalanceCSV << " flux(" << phihwater.boundaryField()[patchi].patch().name() << ")";
            }
        }
        CmassBalanceCSV << " flux(fixedPoint)" << endl;
    }
}
