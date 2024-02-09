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

///////////////////////////////////////////////////////////////////
///////////////////// VELOCITY - FLUXES ///////////////////////////
///////////////////////////////////////////////////////////////////

Info << nl << "Creating fluidPhase" << endl;
word phaseName(transportProperties.optionalSubDict("porousTransport").lookupOrDefault<word>("phaseName",""));
autoPtr<fluidPhase> phase = fluidPhase::New(mesh, transportProperties, phaseName);
volVectorField& Utheta = phase->U();
surfaceScalarField& phi = phase->phi();

//- list that receives event files for active event-based boundary conditions
List<patchEventFile*> patchEventList;
eventFlux::setEventFileRegistry(&patchEventList, "C");

porousMediumModel pmModel(mesh, transportProperties, "");
pmModel.check_eps();
const volScalarField& eps = pmModel.eps();
autoPtr<porousMediumTransportModel> pmTransportModel = porousMediumTransportModel::New(phaseName, pmModel);
multiscalarMixture& composition = pmTransportModel->composition();
List<sourceEventFile*>& tracerSourceEventList = pmTransportModel->sourceEventList();


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
    Info << "===> field S" << phaseName << " computed from theta/eps, min = " << gMin(Saturation.internalField()) << " ; max = " << gMax(Saturation.internalField()) << endl;
}
else
{
    if (Saturation.headerOk())
    {
        theta = Saturation * eps;
        Info << "===> field S" << phaseName << " read, min = " << gMin(Saturation.internalField()) << " ; max = " << gMax(Saturation.internalField()) <<  endl;
    }
    else
    {
        Info << "===> Saturation/Water content file not found (saturated flow)" << endl;
    }
}
