{
    hwater.storePrevIter();

    fvScalarMatrix hwaterEqn
        (
            fvm::div(phi,hwater,"div(phi,hwater)")
            ==
            infiltration
        );
    if (sousrelax < 0.99)
    {
        hwaterEqn.relax(sousrelax);
    }
    else
    {
        hwaterEqn.relax(0.99);
    }
    hwaterEqn.solve();
    hwater.relax();
    //- updating potential
    potential == hwater + z0;
    resPicard = gMax((mag(hwater-hwater.prevIter()))->internalField());
    Info << "Potential min : " << min(potential).value() << " max = " << max(potential).value() << " residual = " << resPicard << endl;

    //- updating flux/velocity field
    phi == (-Mf * fvc::snGrad(potential)) * mesh.magSf();
    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(U.boundaryField()[patchi]))
        {
            phi.boundaryFieldRef()[patchi] = U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

}
