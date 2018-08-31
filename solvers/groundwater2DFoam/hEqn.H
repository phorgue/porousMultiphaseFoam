{

    hwater.storePrevIter();

    fvScalarMatrix hwaterEqn
        (
            eps * fvm::ddt(hwater)
            + fvm::div(phi,hwater,"div(phi,hwater)")
            ==
            infiltration
        );

    hwaterEqn.solve();
 
    dhdTmax = max(mag(fvc::ddt(hwater))).value();

    //- updating potential
    potential == hwater + z0;
    Info << "Potential min : " << min(potential).value() << " max = " << max(potential).value() << " max(dh) = " << dhdTmax*runTime.deltaTValue() << endl;

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
