////////////////////////////////////////////////////////////////
////////////////////// GENERAL PROPERTIES //////////////////////
////////////////////////////////////////////////////////////////

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

Info<< "Reading transportProperties\n" << endl;
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
Info << "}" << endl;

//////////////////////////////////////////////////////////////
////////////////////// PHASE PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////

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

if (theta.headerOk())
{
    Saturation == theta/eps;
    Info << "===> field S" << phaseName << " computed from theta/eps, min = " << min(Saturation).value() << " ; max = " << max(Saturation).value() << nl << endl;
}
else
{
    if (Saturation.headerOk())
    {
        Info << "===> field S" << phaseName << " read, min = " << min(Saturation).value() << " ; max = " << max(Saturation).value() << nl << endl;
    }
    else
    {
        Info << "===> Saturation file not found (saturated flow)" << endl;
    }
}

volVectorField U
(
    IOobject
    (
        "U"+phaseName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Udefault",dimLength/dimTime,vector(0,0,0))
);

Info<< nl << "Reading/Computing field phi" << phaseName << endl;
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
    Info<< "Reading field phi" << phaseName << nl << endl;
}
else
{
    if (U.headerOk())
    {
        Info<< "Computing field phi" << phaseName << " from field U" << phaseName << nl << endl;
    }
    else
    {
        FatalErrorIn("createFields.H") << " neither phi" << phaseName << " or U" << phaseName << " field is present in time directory" << abort(FatalError);
    }
}

//- creation of dispersion model
autoPtr<dispersionModel> DeffModel= dispersionModel::New("DeffModel",transportProperties,U);

Info << "Computing retard coefficient R" << endl;
scalar R =  1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps.value();

////////////////////////////////////////////////////
//////////////////// OUTPUT CSV ////////////////////
////////////////////////////////////////////////////

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
            massBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
        }
    }
    massBalanceCSV << endl;
}
