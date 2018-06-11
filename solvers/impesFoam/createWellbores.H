Info<< "Reading wellboreProperties\n" << endl;

dictionary wellboreProperties(transportProperties.subDict("wellboreProperties"));

// Injection Wellbore
dimensionedScalar qInjection
(
    wellboreProperties.lookupOrDefault<dimensionedScalar>
    (
        "qInjection",
        dimensionedScalar("qInjection",dimVol/dimTime,0.)
    )
);    

// Extraction Wellbore
dimensionedScalar qExtraction
(
    wellboreProperties.lookupOrDefault<dimensionedScalar>
    (
        "qExtraction",
        dimensionedScalar("qExtraction",dimVol/dimTime,0.)
    )
);

Info << " injection flow (m3/s) = " << qInjection << endl;
Info << " extraction flow (m3/s) = " << qExtraction << endl;

Info<< "Reading field Winj (injection wellbore)\n" << endl;
volScalarField Winj
(
    IOobject
    (
        "Winj",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimensionSet(0,0,0,0,0),0)
);

dimensionedScalar Vinj(sum(Winj*mesh.V()));
dimensionedScalar SrcInj("SrcInj",dimensionSet(0,0,-1,0,0),0.);
SrcInj = qInjection/(Vinj+dimensionedScalar("",dimVol,SMALL));

Info<< "Reading field Wext (extraction wellbore)\n" << endl;
volScalarField Wext
(
    IOobject
    (
        "Wext",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimensionSet(0,0,0,0,0),0)
);

dimensionedScalar Vext(sum(Wext*mesh.V()));
dimensionedScalar SrcExt("SrcExt",dimensionSet(0,0,-1,0,0),0.);
SrcExt =  qExtraction/(Vinj+dimensionedScalar("",dimVol,SMALL));
