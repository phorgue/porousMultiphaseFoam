{
        //- update dispersion tensor coefficient
        DeffModel->correct();
        
        fvScalarMatrix CEqn
            (
                eps * R *fvm::ddt(C)
                + fvm::div(phi, C, "div(phi,C)")
                - fvm::laplacian(DeffModel->Deff(), C, "laplacian(Deff,C)")
                ==
                - eps * R * fvm::Sp(lambda,C)
            );

        CEqn.solve();
        dCdTmax = max(mag(fvc::ddt(C))).value();
        Info << "Concentration Min(C) = " << min(C).value() << " Max(C) = " << max(C).value() << " dCdTmax = " << dCdTmax << endl;

}
