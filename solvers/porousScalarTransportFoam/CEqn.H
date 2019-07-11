{
    //- update dispersion tensor coefficient
    DeffModel->correct(C,U,Saturation,eps);

    volTensorField Deff = eps*Saturation*DeffModel->Deff();

    fvScalarMatrix CEqn
        (
            eps * R * Saturation * fvm::ddt(C)
            + fvm::div(phi, C, "div(phi,C)")
            - fvm::laplacian(Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm
            - eps * R * Saturation * fvm::Sp(lambda,C)
        );

    CEqn.solve();

    dCdTmax = max(mag(fvc::ddt(C))).value();

    Info << "Concentration Min(C) = " << min(C).value() << " Max(C) = " << max(C).value() << " mass(C) = " << fvc::domainIntegrate(R*C*Saturation*eps).value() << endl;

}
