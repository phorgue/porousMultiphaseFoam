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
    transportProperties.lookup("eps")
);

Info<< "Reading porous transport properties" << endl;
dictionary porousTransport(transportProperties.subDict("porousTransport"));
word phaseName(porousTransport.lookupOrDefault<word>("phaseName","notDefined"));
Info << nl << "porousTransport parameters" << nl << "{" << endl;
Info << "    " << "phaseName " << " : " << phaseName << endl;
Info << "} \n" << endl;


Info<< "Reading field U" << phaseName << nl << endl;
volVectorField U
(
    IOobject
    (
        "U"+phaseName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field phi" << phaseName << nl << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi"+phaseName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- creation of dispersion model
autoPtr<dispersionModel> DeffModel= dispersionModel::New("DeffModel",transportProperties,U);
