{
    potential.storePrevIter();

    fvScalarMatrix potentialEqn
        (
            eps * fvm::ddt(potential)
            - fvm::laplacian(transmissivity,potential,"laplacian(transmissivity,potential)")
            ==
            - waterSourceTerm * zScale
            - infiltration
        );

    #include "updateForcing.H"

    potentialEqn.solve();

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
        Info << "Number of dry cells = " << dryCellIDList.size() << ", water added = " << waterAdded << " m3, cumulative water added = " << cumulativeWaterAdded << " m3 (" << cumulativeWaterAdded*zScale/sum(mesh.V()).value() << " m)" << endl;
    }

    //- compute derivatives
    dtManager.updateDerivatives();

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

    Info << "Potential min: " << gMin(potential.internalField()) << ", max = " << gMax(potential.internalField()) << ", delta(potential) = " << dtManager.dVmax()*runTime.deltaTValue() << endl;

}
