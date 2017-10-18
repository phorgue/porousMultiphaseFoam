/////////////////////////////////////////////////////////////////////////////
//////////////////////////// MESH.delta()  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// - compute fvMesh.delta() explicitly to be compatible with Foam-extend
surfaceVectorField delta
(
    IOobject
    (
        "delta",
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimLength
);
const volVectorField& C = mesh.C();
const unallocLabelList& owner = mesh.owner();
const unallocLabelList& neighbour = mesh.neighbour();


forAll(owner, facei)
{
    delta[facei] = C[neighbour[facei]] - C[owner[facei]];
}
forAll(delta.boundaryField(), patchi)
{
    delta.boundaryField()[patchi] = mesh.boundary()[patchi].delta();
}
