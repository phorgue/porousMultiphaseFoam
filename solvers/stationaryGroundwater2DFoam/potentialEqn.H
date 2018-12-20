{
    potential.storePrevIter();

    //- computing seep term if active
    if (forcingPotential == "seepage")
    {
        forAll(mesh.C(),celli)
        {
            seepTerm[celli] = max(0,(potential[celli]-potentialMNT[celli]));
        }
        Info << "seep term : min = " << min(seepTerm.internalField()).value() << ", max = " << max(seepTerm.internalField()).value() << endl;;
    }

    fvScalarMatrix potentialEqn
        (
            - fvm::laplacian(transmissivity,potential,"laplacian(transmissivity,potential)")
            ==
            infiltration - seepTerm
        );

    //- Fix values if forcingPotential "fixedValue" or "fixedMNT"
    if ((forcingPotential.substr(0,5) == "fixed") && (fixedPotentialIDList.size() > 0))
    {
        potentialEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);
    }

    potentialEqn.solve();
    potential.relax();

    Info << "Potential min : " << min(potential).value() << ", max = " << max(potential).value() << ", diff = " << max(mag(potential-potential.oldTime())).value() << endl;

    //- updating flow properties
    hwater == potential - z0;
    if (min(hwater).value() <= 0) FatalErrorIn("potentialEqn.H") << " Computed hwater fields has negative values, try to reduce relaxation" << exit(FatalError);

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
