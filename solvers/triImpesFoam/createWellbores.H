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

// WATER INJECTION 
//- injection Wellbore
dimensionedScalar qInjectionw(wellboreProperties.lookup("qInjectionw"));    

Info<< "Reading field Winjw (injection wellbore of water)\n" << endl;
volScalarField Winjw
(
    IOobject
    (
        "Winjw",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimensionSet(0,0,0,0,0),0)
);

dimensionedScalar Vinjw(sum(Winjw*mesh.V()));
dimensionedScalar SrcInjw("SrcInjw",dimensionSet(0,0,-1,0,0),0.);
SrcInjw = qInjectionw/(Vinjw+dimensionedScalar("",dimVol,SMALL));

// Extraction Wellbore
dimensionedScalar qExtractionw
( 
    wellboreProperties.lookup("qExtractionw")  
);    

Info<< "Reading field Wextw (extraction wellbore of water)\n" << endl;
volScalarField Wextw
(
    IOobject
    (
        "Wextw",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimensionSet(0,0,0,0,0),0)
);

dimensionedScalar Vextw(sum(Wextw*mesh.V()));
dimensionedScalar SrcExtw("SrcExtw",dimensionSet(0,0,-1,0,0),0.);
SrcExtw = qExtractionw/(Vinjw+dimensionedScalar("",dimVol,SMALL));

// OIL INJECTION 
//- injection Wellbore
dimensionedScalar qInjectiono(wellboreProperties.lookup("qInjectiono"));    

Info<< "Reading field Winjo (injection wellbore of oil)\n" << endl;
volScalarField Winjo
(
    IOobject
    (
        "Winjo",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimensionSet(0,0,0,0,0),0)
);

dimensionedScalar Vinjo(sum(Winjo*mesh.V()));
dimensionedScalar SrcInjo("SrcInjo",dimensionSet(0,0,-1,0,0),0.);
SrcInjo = qInjectiono/(Vinjo+dimensionedScalar("",dimVol,SMALL));

// Extraction Wellbore
dimensionedScalar qExtractiono
( 
    wellboreProperties.lookup("qExtractiono")  
);    

Info<< "Reading field Wexto (extraction wellbore of oil)\n" << endl;
volScalarField Wexto
(
    IOobject
    (
        "Wexto",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimensionSet(0,0,0,0,0),0)
);

dimensionedScalar Vexto(sum(Wexto*mesh.V()));
dimensionedScalar SrcExto("SrcExto",dimensionSet(0,0,-1,0,0),0.);
SrcExto = qExtractiono/(Vinjo+dimensionedScalar("",dimVol,SMALL));
