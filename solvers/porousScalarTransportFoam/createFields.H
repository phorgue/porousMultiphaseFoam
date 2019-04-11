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
        IOobject::MUST_READ_IF_MODIFIED,
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

word Uname="U"+phaseName;
if (theta.headerOk())
{
    Saturation == theta/eps;
    Uname = "U";
    Info << "===> field S" << phaseName << " computed from theta/eps, min = " << min(Saturation).value() << " ; max = " << max(Saturation).value() << endl;
}
else
{
    if (Saturation.headerOk())
    {
        Info << "===> field S" << phaseName << " read, min = " << min(Saturation).value() << " ; max = " << max(Saturation).value() <<  endl;
    }
    else
    {
        Info << "===> Saturation file not found (saturated flow)" << endl;
    }
}

Info<< nl << "Reading field U" << phaseName << endl;
volVectorField U
(
    IOobject
    (
        Uname,
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
    Info<< nl << "Computing field phi" << phaseName << " from field U" << phaseName << endl;
}

//- creation of dispersion model
autoPtr<dispersionModel> DeffModel= dispersionModel::New("DeffModel",transportProperties,mesh);

Info << "Computing retard coefficient R = ";
scalar R =  1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps.value();
Info << R << endl;

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
    dimensionedScalar("",dimMass/(dimVol*dimTime),0)
);

//- Looking for 2nd order time-scheme for source term
bool isBackwardScheme = false;
if (word(mesh.ddtScheme("source")) == "backward")
{
    isBackwardScheme = true;
}
