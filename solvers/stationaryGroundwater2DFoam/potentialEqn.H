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

    //- updating flow properties
    hwater = potential - z0;

    //- Checking for dry cells
    if (min(hwater).value() <= hwaterMin.value())
    {
        label ndryCells = 0;
        forAll(hwater, celli)
        {
            if (hwater[celli] <= hwaterMin.value())
            {
                ndryCells++;
                hwater[celli] = hwaterMin.value();
            }
        }
        Warning() << "Number of dry cells = " << ndryCells << endl;
    }

    Info << "Potential min : " << min(potential).value() << ", max = " << max(potential).value() << ", diff = " << max(mag(potential-potential.oldTime())).value() << endl;

    transmissivity = Mf*fvc::interpolate(hwater);
    phi = (-Mf * fvc::snGrad(potential)) * mesh.magSf();
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
