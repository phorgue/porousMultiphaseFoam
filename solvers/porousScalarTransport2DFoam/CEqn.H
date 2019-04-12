{
    //- update dispersion tensor coefficient
    DeffModel->correct(C,U,eps);

    volTensorField Deff = eps*hwater*DeffModel->Deff();

    fvScalarMatrix CEqn
        (
            eps * R * hwater * fvm::ddt(C)
            + fvm::div(phiWater, C, "div(phi,C)")
            - fvm::laplacian(Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm * zScale
            - eps * R * hwater * fvm::Sp(lambda,C)
            - dimensionedScalar("",dimLength,zScale) * fvm::Sp(seepageTerm,C)
        );

    CEqn.solve();

    Info << "Concentration Min(C) = " << min(C).value() << " Max(C) = " << max(C).value()  << " mass(C) = " << fvc::domainIntegrate(R*C*hwater*eps).value()/zScale << endl;

}
