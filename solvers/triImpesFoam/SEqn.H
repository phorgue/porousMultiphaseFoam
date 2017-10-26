{

    Ua.correctBoundaryConditions();
    Uo.correctBoundaryConditions();
    Uw.correctBoundaryConditions();   
    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(Ua.boundaryField()[patchi]))
        {
            phia.boundaryFieldRef()[patchi] = Ua.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
        if (isA< fixedValueFvPatchField<vector> >(Uo.boundaryField()[patchi]))
        {
            phio.boundaryFieldRef()[patchi] = Uo.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
        if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]))
        {
            phiw.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }

    fvScalarMatrix SoEqn
        (
            eps*fvm::ddt(So) + fvc::div(phio) 
            // wellbore part
            - (SrcInjo*Winjo-Fo*SrcExto*Wexto)
        );

    SoEqn.solve();

    Info << "Saturation oil " << " Min(So) = " << gMin(So) << " Max(So) = " << gMax(So) << " Volume = " << fvc::domainIntegrate(eps*So).value() <<  endl;

    fvScalarMatrix SwEqn
        (
            eps*fvm::ddt(Sw) + fvc::div(phiw) 
            // wellbore part
            - (SrcInjw*Winjw-Fw*SrcExtw*Wextw)
        );

    SwEqn.solve();

    Info << "Saturation water " << " Min(Sw) = " << gMin(Sw) << " Max(Sw) = " << gMax(Sw) << " Volume = " << fvc::domainIntegrate(eps*Sw).value() <<  endl;

    volScalarField Stotal = Sw+So;
    Info << "Saturation total " << " Min(Stotal) = " << gMin(Stotal) << " Max(Stotal) = " << gMax(Stotal) << endl;
        
}
