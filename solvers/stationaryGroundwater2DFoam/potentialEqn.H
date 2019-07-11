{
    potential.storePrevIter();

    fvScalarMatrix potentialEqn
        (
            - fvm::laplacian(transmissivity,potential,"laplacian(transmissivity,potential)")
            ==
            - infiltration
        );

    #include "updateForcing.H"

    potentialEqn.solve();
    potential.relax();
    Info << "Potential min : " << min(potential).value() << ", max = " << max(potential).value() << ", diff = " << max(mag(potential-potential.oldTime())).value() << endl;

    //- updating flow properties
    hwater == potential - z0;
    if (min(hwater).value() <= 0) FatalErrorIn("potentialEqn.H") << " Computed hwater fields has negative values" << exit(FatalError);

    transmissivity = Mf*fvc::interpolate(hwater);
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
