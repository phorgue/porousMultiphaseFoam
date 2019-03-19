{
    //- update dispersion tensor coefficient
    DeffModel->correct(C,U,eps);

    fvScalarMatrix CEqn
        (
            eps * R * hwater * fvm::ddt(C)
            + fvm::div(phi, C, "div(phi,C)")
            - eps * hwater * fvm::laplacian(DeffModel->Deff(), C, "laplacian(Deff,C)")
            ==
            - R * sourceTerm * zScale
            - eps * R * hwater * fvm::Sp(lambda,C)
        );

    CEqn.solve();

    Info << "Concentration Min(C) = " << min(C).value() << " Max(C) = " << max(C).value()  << " mass(C) = " << fvc::domainIntegrate(C*hwater*eps).value()/zScale << endl;

}
