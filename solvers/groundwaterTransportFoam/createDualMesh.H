autoPtr<dynamicFvMesh> meshPtrFluid(dynamicFvMesh::New(args, runTime));
dynamicFvMesh& mesh = meshPtrFluid();
mesh.init(true);   // initialise all (lower levels and current)

dynamicFvMesh* meshPtrTransport(&mesh);
if (dualMesh)
{
    autoPtr<dynamicFvMesh> meshAutoPtrTransport = dynamicFvMesh::New
            (
                    IOobject
                            (
                                    "transport",
                                    runTime.timeName(),
                                    runTime,
                                    IOobject::MUST_READ
                            )
            );
    meshPtrTransport = meshAutoPtrTransport.release();
}
dynamicFvMesh& meshT = *meshPtrTransport;
meshT.init(true);

//const labelIOList& cell_level =  refCast<dynamicRefineFvMesh>(meshT).meshCutter().cellLevel();
//const refinementHistory& ref_hist = refCast<dynamicRefineFvMesh>(meshT).meshCutter().history();

//- To protect boundary cell from refining
if (meshT.dynamic()) {
    if (isA<dynamicRefineFvMesh>(meshT)) {
        DynamicList<label> boundary_protected_cells(0);
        const polyBoundaryMesh& patches = meshT.boundaryMesh();
        forAll(patches, patchi) {
            if (not(isA<processorPolyPatch>(patches[patchi])) &&
                not(isA<symmetryPlanePolyPatch>(patches[patchi])) ){
                boundary_protected_cells.append(meshT.boundary()[patchi].faceCells());
            }
        }
        refCast<dynamicRefineFvMesh>(meshT).protectedCell() = bitSet(boundary_protected_cells);
        }
    DynamicList<label> boundary_protected_cells(0);
}
