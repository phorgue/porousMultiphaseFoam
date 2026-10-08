{

    Ua.correctBoundaryConditions();
    Ub.correctBoundaryConditions();
    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(Ua.boundaryField()[patchi]))
        {
            phia.boundaryFieldRef()[patchi] = Ua.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
        if (isA< fixedValueFvPatchField<vector> >(Ub.boundaryField()[patchi]))
        {
            phib.boundaryFieldRef()[patchi] = Ub.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
  
    fvScalarMatrix SbEqn
        (
            // eps*fvm::ddt(Sb) + fvc::div(phib)
            eps*fvm::ddt(Sb) + fvm::div(phib, Sb) 
            ==
            // event source terms
            - sourceTerm
        );

    SbEqn.solve();
    Sb.correctBoundaryConditions();

    Info << "Saturation b: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

}
