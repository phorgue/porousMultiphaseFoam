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
        dryCellIDList.clear();
        forAll(hwater, celli)
        {
            if (hwater[celli] <= hwaterMin.value())
            {
                dryCellIDList.append(celli);
                hwater[celli] = hwaterMin.value();
            }
        }
        Info << "Number of dry cells = " << dryCellIDList.size() << endl;
    }

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
    cellFlux = fvc::div(phi*fvc::interpolate(hwater)) + infiltration;

    Info << "Potential min : " << min(potential).value() << ", max = " << max(potential).value() << ", diff = " << max(mag(potential-potential.oldTime())).value() << endl;

}
