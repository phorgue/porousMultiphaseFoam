{
    potential.storePrevIter();

    fvScalarMatrix potentialEqn
        (
            - fvm::laplacian(transmissivity,potential,"laplacian(transmissivity,potential)")
            ==
            - waterSourceTerm * zScale
            - infiltration
        );

    if (!steady) potentialEqn += eps * fvm::ddt(potential);

    #include "updateForcing.H"

    potentialEqn.solve();

    if (steady) potential.relax();

    //- updating hwater and fluxes
    hwater = potential - z0;

    //- Checking for dry cells
    if (gMin(hwater.internalField()) <= hwaterMin.value())
    {
        scalar waterAdded = 0;
        dryCellIDList.clear();
        forAll(hwater,celli)
        {
            if (hwater[celli] <= hwaterMin.value())
            {
                dryCellIDList.append(celli);
                waterAdded += (hwaterMin.value()-hwater[celli])*mesh.V()[celli]/zScale;
                hwater[celli] = hwaterMin.value();
            }
        }
        cumulativeWaterAdded += waterAdded;
        Info << "Number of dry cells = " << dryCellIDList.size();
        if (steady)
        {
            Info << ", water added = " << waterAdded << " m3, cumulative water added = " << cumulativeWaterAdded << " m3 (" << cumulativeWaterAdded*zScale/gSum(mesh.V()) << " m)";
        }
        Info << endl;
    }

    //- compute derivatives
    if (!steady) dtManager.updateDerivatives();

    //- updating flow properties
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
    cellFlux = fvc::div(phi*fvc::interpolate(hwater)) + infiltration + zScale * waterSourceTerm;

    Info << "Potential min: " << gMin(potential.internalField()) << ", max = " << gMax(potential.internalField()) << ", delta(potential) = " << gMax(mag(potential.internalField()-potential.oldTime().internalField())()) << endl;
}
