{
        //- update dispersion tensor coefficient
        DeffModel->correct();
        
        fvScalarMatrix CEqn
            (
                eps*fvm::ddt(C)
                + fvm::div(phi, C)
                - fvm::laplacian(DeffModel->Deff(), C, "laplacian(Deff,C)")
            );

        CEqn.solve();
        dCdTmax = max(mag(fvc::ddt(C))).value();
        Info << "Min(C) = " << min(C).value() << " Max(C) = " << max(C).value() << " dCdTmax = " << dCdTmax << endl;

}
