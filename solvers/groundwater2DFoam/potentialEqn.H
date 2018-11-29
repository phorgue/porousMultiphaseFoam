{
    potential.storePrevIter();

    fvScalarMatrix potentialEqn
        (
            eps * fvm::ddt(potential)
            - fvm::laplacian(transmissivity,potential,"laplacian(transmissivity,potential)")
            ==
            infiltration
        );

    if(fixedPotentialIDList.size() > 0) potentialEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);

    potentialEqn.solve();

    //- updating hwater
    hwater == potential - z0;
    if (min(hwater).value() <= 0) FatalErrorIn("potentialEqn.H") << " Computed hwater fields has negative values, try to reduce relaxation" << exit(FatalError);
    dhdTmax = (max(mag(fvc::ddt(hwater)))/max(hwater)).value();
    Info << "Potential min : " << min(potential).value() << ", max = " << max(potential).value() << ", max(dh) = " << dhdTmax*runTime.deltaTValue() << endl;

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
