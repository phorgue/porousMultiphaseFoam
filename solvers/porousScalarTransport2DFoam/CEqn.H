{
    //- update dispersion tensor coefficient
    DeffModel->correct(C,U);

    fvScalarMatrix CEqn
        (
            eps * R * hwater * fvm::ddt(C)
            + dimensionedScalar("",dimLength,1) * fvm::div(phi, C, "div(phi,C)")
            - eps * hwater * fvm::laplacian(DeffModel->Deff(), C, "laplacian(Deff,C)")
            ==
            R * sourceTerm * dimensionedScalar("",dimLength,1)
            - eps * R * hwater * fvm::Sp(lambda,C)
        );

    CEqn.solve();

    dCdTmax = max(mag(fvc::ddt(C))).value();

    Info << "Concentration Min(C) = " << min(C).value() << " Max(C) = " << max(C).value() << " dCdTmax = " << dCdTmax*runTime.deltaTValue() << " mass(C) = " << fvc::domainIntegrate(C*hwater*eps).value() << endl;

}
