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
    hwater = potential - z0;

    //- compute derivatives
    dhwaterdTmax = max(mag(fvc::ddt(hwater))).value();
    if (timeScheme == "Euler")
    {
        volScalarField dhwater2dT2(d2dt2Operator.fvcD2dt2(hwater));
        dhwater2dT2max = 0;
        forAll(dhwater2dT2, celli)
        {
            if(mag(dhwater2dT2[celli]) > dhwater2dT2max)
            {
                hwatermax = hwater[celli];
                dhwater2dT2max = mag(dhwater2dT2[celli]);
            }
        }
    }
    else
    {
        volScalarField dhwater3dT3(d3dt3Operator.d3dt3(hwater));
        dhwater3dT3max = 0;
        forAll(dhwater3dT3, celli)
        {
            if(mag(dhwater3dT3[celli]) > dhwater3dT3max)
            {
                hwatermax = hwater[celli];
                dhwater3dT3max = mag(dhwater3dT3[celli]);
            }
        }
    }

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
        Warning() << "number of dry cells = " << ndryCells << ", water added = " << waterAdded << " m3, cumulative water added = " << cumulativeWaterAdded << " m3 (" << cumulativeWaterAdded*zScale/sum(mesh.V()).value() << " m)" << endl;
     }

    Info << "Potential min : " << min(potential).value() << ", max = " << max(potential).value() << ", max(dhwater) = " << dhwaterdTmax*runTime.deltaTValue() << endl;

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

}
