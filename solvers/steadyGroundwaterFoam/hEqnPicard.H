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

    #include "updateForcing.H"

    hEqn.solve();

    deltah = mag(h-h.prevIter());
}
