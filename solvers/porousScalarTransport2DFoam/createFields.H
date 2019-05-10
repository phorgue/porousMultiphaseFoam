////////////////////////////////////////////////////////////////
////////////////////// GENERAL PROPERTIES //////////////////////
////////////////////////////////////////////////////////////////

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

Info<< "Reading porosity eps\n" << endl;
dimensionedScalar eps
(
    transportProperties.lookupOrDefault("eps",dimensionedScalar("eps_default",dimless,1))
);

Info<< "Reading porous transport properties" << endl;
dictionary porousTransport(transportProperties.subDict("porousTransport"));
word phaseName(porousTransport.lookupOrDefault<word>("phaseName",""));
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

//////////////////////////////////////////////////////////////
////////////////////// PHASE PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////

Info << "Reading hwater field if present...";
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
    Info << nl << "min(hwater) = " << min(hwater).value() << " ; max(hwater) = " << max(hwater).value() << endl;
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
    Info << "min(z0) = " << min(z0).value() << " ; max(z0) = " << max(z0).value() << endl;

    hwater == potential - z0;

    Info << nl << "Computed hwater : min(hwater) = " << min(hwater).value() << " ; max(hwater) = " << max(hwater).value() << endl;
    hwater.write();
}

if (min(hwater).value() < 0)
{
    FatalErrorIn("createFields") << " hwater has negative values" << abort(FatalError);
}

Info << nl << "Reading field U" << phaseName << endl;
volVectorField U
(
    IOobject
    (
        "U"+phaseName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

surfaceScalarField phi
(
    IOobject
    (
        "phi"+phaseName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);

if (phi.headerOk())
{
    Info << nl << "Reading field phi" << phaseName << endl;
}
else
{
    Info<< "Computing field phi" << phaseName << " from field U" << phaseName << nl << endl;
    phi.write();
}

surfaceScalarField phiWater = phi * fvc::interpolate(hwater);

//- creation of dispersion model
autoPtr<dispersionModel> DeffModel= dispersionModel::New("DeffModel",transportProperties, mesh);
volTensorField Deff =  eps*hwater*DeffModel->Deff();

Info << "Computing retard coefficient R = ";
scalar R =  1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps.value();
Info << R << endl;

/////////////////////////////////////////////////////////
//////////////////// OUTLET CSV FILE ////////////////////
/////////////////////////////////////////////////////////

bool CSVactive=runTime.controlDict().lookupOrDefault<bool>("CSVactive",true);
word CSVoutput=runTime.controlDict().lookupOrDefault<word>("CSVoutput","CmassBalance.csv");
OFstream massBalanceCSV(CSVoutput);
if (CSVactive)
{
    massBalanceCSV << "#Time TotalMass(kg)";
    forAll(mesh.boundaryMesh(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            massBalanceCSV << " flux(" << phiWater.boundaryField()[patchi].patch().name() << ")";
        }
    }
    massBalanceCSV << " flux(fixedPoint)" << endl;
}

//- Source term
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

//- Looking for 2nd order time-scheme for source term
bool isBackwardScheme = false;
if (word(mesh.ddtScheme("source")) == "backward")
{
    isBackwardScheme = true;
}
