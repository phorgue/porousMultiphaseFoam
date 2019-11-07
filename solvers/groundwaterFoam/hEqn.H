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
            //- accumulation terms
            (Ss*pcModel->Se() + pcModel->Ch()) * fvm::ddt(h)
            //-mass conservative terms
            + massConservativeTerms * (
                pcModel->Ch()*(h.oldTime()-h.prevIter())
                + ( theta - theta.oldTime()))
            /runTime.deltaT()
            //- transport terms
            - fvm::laplacian(Mf,h)
            + fvc::div(phiG)
            ==
            - sourceTerm
        );

    #include "updateForcing.H"

    hEqn.solve();

    resPicard = gMax((mag(h-h.prevIter()))->internalField());
    Info << "Picard iteration : " << iterPicard << " residual = " << resPicard << endl; 

    phi = phiG-(Mf*fvc::snGrad(h))*mesh.magSf();
    Utheta = fvc::reconstruct(phi);

}
