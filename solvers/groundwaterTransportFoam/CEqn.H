{
        //- update dispersion tensor coefficient
        DeffModel->correct();
        
        fvScalarMatrix CEqn
            (
                R * theta * fvm::ddt(C)
                + R * C * fvc::ddt(theta)
                + fvm::div(phi, C, "div(phi,C)")
                - fvm::laplacian(DeffModel->Deff(), C, "laplacian(Deff,C)")
                ==
                R * (SrcInj*Cin)*Winj - (SrcExt*C)*Wext
                - R * theta * fvm::Sp(lambda,C)
            );

        CEqn.solve();

        dCdTmax = max(mag(fvc::ddt(C))).value();
        dC = dCdTmax*runTime.deltaTValue();

        Info << "Concentration Min(C) = " << min(C).value() << " Max(C) = " << max(C).value() << " dCdTmax = " << dCdTmax << " mass(C) = " << fvc::domainIntegrate(C*theta).value() << endl;

}
