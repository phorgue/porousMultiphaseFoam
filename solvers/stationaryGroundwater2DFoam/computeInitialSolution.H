{
    //- Compute initial solution
    solve(fvm::laplacian(potential));
    Info << "Potential min : " << min(potential).value() << " max = " << max(potential).value() << endl;
    
    //- Updating hwater and phi fields
    hwater == potential - z0;

    phi == (-Mf * fvc::snGrad(potential)) * mesh.magSf();
    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(U.boundaryField()[patchi]))
        {
            phi.boundaryFieldRef()[patchi] = U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
}
