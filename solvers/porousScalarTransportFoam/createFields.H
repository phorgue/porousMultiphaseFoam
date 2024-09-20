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
autoPtr<porousMediumTransportModel> pmTransportModel = porousMediumTransportModel::New(phaseName, mesh, transportProperties);
multiscalarMixture& composition = pmTransportModel->composition();
composition.initRetardCoef(eps);
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
    Info << "===> field theta read , min = " << gMin(theta.internalField()) << " ; max = " << gMax(theta.internalField()) <<  endl;
}
else
{
    if (Saturation.headerOk())
    {
theta = Saturation * eps;
Info << "===> field theta" << phaseName << " computed, min = " << gMin(theta.internalField()) << " ; max = " << gMax(theta.internalField()) <<  endl;
    }
    else
    {
        Info << "===> Saturation/Water content file not found (saturated flow)" << endl;
    }
}

//- To protect boundary cell from refining
if (isA<dynamicRefineFvMesh>(mesh)) {
    DynamicList<label> boundary_protected_cells(0);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchi) {
        if (not(isA<processorPolyPatch>(patches[patchi])) &&
                not(isA<symmetryPlanePolyPatch>(patches[patchi])) ){
            boundary_protected_cells.append(mesh.boundary()[patchi].faceCells());
        }
    }
    refCast<dynamicRefineFvMesh>(mesh).protectedCell() = bitSet(boundary_protected_cells);
}

