Info<< "Reading wellboreProperties\n" << endl;

IOdictionary wellboreProperties
(
    IOobject
    (
        "wellboreProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar activateWellbores(wellboreProperties.lookupOrDefault<scalar>("activateWellbores",0.)); 

// Injection Wellbore
dimensionedScalar qInjection(wellboreProperties.lookup("qInjection"));    

Info<< "Reading field Winj (injection wellbore)\n" << endl;
volScalarField Winj
(
    IOobject
    (
        "Winj",
        runTime.timeName(),
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

// Extraction Wellbore
dimensionedScalar qExtraction
( 
    wellboreProperties.lookup("qExtraction")  
);    

Info<< "Reading field Wext (extraction wellbore)\n" << endl;
volScalarField Wext
(
    IOobject
    (
        "Wext",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimensionSet(0,0,0,0,0),0)
);

dimensionedScalar Vext(sum(Wext*mesh.V()));
dimensionedScalar SrcExt("SrcExt",dimensionSet(0,0,-1,0,0),0.);
SrcExt = qExtraction/(Vinj+dimensionedScalar("",dimVol,SMALL));
