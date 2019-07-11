{

    Utheta.correctBoundaryConditions();
    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(Utheta.boundaryField()[patchi]))
        {
            phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
    
    h.storePrevIter();

    fvScalarMatrix hEqn
        (
            - fvm::laplacian(Mf,h)
            + fvc::div(phiG)
            ==
            - sourceTerm
        );

    hEqn.solve();
    h.relax();

    phi == phiG-(Mf*fvc::snGrad(h))*mesh.magSf();
    Utheta = fvc::reconstruct(phi);

}
