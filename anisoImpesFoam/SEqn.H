{

    fvScalarMatrix SbEqn
        (
            eps*fvm::ddt(Sb) + fvc::div(phib) 
            // wellbore part
            - (SrcInj*Winj-Fb*SrcExt*Wext)*activateWellbores
        );

    SbEqn.solve();

    Info << "Saturation b " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) <<  endl;
}
