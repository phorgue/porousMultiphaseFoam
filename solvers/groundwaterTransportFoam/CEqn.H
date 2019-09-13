{
    //- update dispersion tensor coefficient
    DeffModel->correct(C,Utheta,theta);

    volTensorField Deff = theta*DeffModel->Deff();

    fvScalarMatrix CEqn
        (
            R * fvm::ddt(theta,C)
            + fvm::div(phi, C, "div(phi,C)")
            - fvm::laplacian(Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm_tracer
            - R * theta * fvm::Sp(lambda,C)
        );

    CEqn.solve();

    dCdTmax = max(mag(fvc::ddt(C))).value();

    Info << "Concentration Min(C) = " << min(C).value() << " Max(C) = " << max(C).value() << " dCdTmax = " << dCdTmax << " mass(C) = " << fvc::domainIntegrate(C*theta).value() << endl;

}
