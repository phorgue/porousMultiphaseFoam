{
    potential.storePrevIter();

    fvScalarMatrix potentialEqn
        (
            eps * fvm::ddt(potential)
            - fvm::laplacian(transmissivity,potential,"laplacian(transmissivity,potential)")
            ==
            - infiltration
        );

    #include "updateForcing.H"

    potentialEqn.solve();

    //- updating hwater and fluxes
    hwater == potential - z0;

    //- Checking for dry cells
    if (min(hwater).value() <= hwaterMin.value())
    {
        scalar waterAdded = 0;
        label ndryCells = 0;
        forAll(hwater,celli)
        {
            if (hwater[celli] <= hwaterMin.value())
            {
                ndryCells++;
                waterAdded += (hwaterMin.value()-hwater[celli])*mesh.V()[celli]/zScale;
                hwater[celli] = hwaterMin.value();
            }
        }
        cumulativeWaterAdded += waterAdded;
        Info << "*** WARNING *** : number of dry cells = " << ndryCells << ", water added = " << waterAdded << " m3, cumulative water added = " << cumulativeWaterAdded << " m3 (" << cumulativeWaterAdded*zScale/sum(mesh.V()).value() << " m)" << endl;
     }

    Info << "Potential min : " << min(potential).value() << ", max = " << max(potential).value() << ", max(dh) = " << dhdTmax*runTime.deltaTValue() << endl;

    dhdTmax = (max(mag(fvc::ddt(hwater)))/max(hwater)).value();

    //- updating flow properties
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
