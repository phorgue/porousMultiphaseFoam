{
        //- update dispersion tensor coefficient
        DeffModel->correct(C);

        fvScalarMatrix CEqn
            (
                eps * R * hwater * fvm::ddt(C)
                + fvm::div(phi, C, "div(phi,C)")
                - fvm::laplacian(hwater*DeffModel->Deff(), C, "laplacian(Deff,C)")
                ==
                R * (SrcInj*Cin)*Winj - (SrcExt*C)*Wext
                - eps * R * hwater * fvm::Sp(lambda,C)
            );

        CEqn.solve();

        dCdTmax = max(mag(fvc::ddt(C))).value();

        Info << "Concentration Min(C) = " << min(C).value() << " Max(C) = " << max(C).value() << " dCdTmax = " << dCdTmax*runTime.deltaTValue() << " mass(C) = " << fvc::domainIntegrate(C*hwater*eps).value() << endl;

}
