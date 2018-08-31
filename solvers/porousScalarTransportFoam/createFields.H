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
word phaseName(porousTransport.lookupOrDefault<word>("phaseName","notDefined"));
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

Info << "Reading Saturation Field (if present)... ";
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
if (Saturation.headerOk())
{
    Info << "field S" << phaseName << " read, min = " << min(Saturation).value() << " ; max = " << max(Saturation).value() << nl << endl;
}
else
{
    Info << "Saturation file not found (saturated flow)" << nl << endl;
}

Info<< "Reading field U" << phaseName << nl << endl;
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

Info<< "Reading/Computing field phi" << phaseName << nl << endl;
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

phi.write(); // write if computed
phi.writeOpt()=IOobject::NO_WRITE;

//- creation of dispersion model
autoPtr<dispersionModel> DeffModel= dispersionModel::New("DeffModel",transportProperties,U);

Info << "Computing retard coefficient R" << endl;
scalar R =  1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps.value();

//////////////////////////////////////////////////////////////
//////////////////// INJECTION PROPERTIES ////////////////////
//////////////////////////////////////////////////////////////

#include "createWellbores.H"
dimensionedScalar Cin(wellboreProperties.lookup("Cin"));
