{

    Ua.correctBoundaryConditions();
    Ub.correctBoundaryConditions();
    forAll(phib.boundaryField(),patchi)
    {
        phia.boundaryFieldRef()[patchi] = Ua.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        phib.boundaryFieldRef()[patchi] = Ub.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    }

    fvScalarMatrix SbEqn
        (
            eps*fvm::ddt(Sb) + fvc::div(phib) 
            // wellbore part
            - (SrcInj*Winj-Fb*SrcExt*Wext)*activateWellbores
        );

    SbEqn.solve();

    Info << "Saturation b " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) <<  endl;

}
