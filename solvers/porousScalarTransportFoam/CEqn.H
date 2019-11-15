forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];
    
    auto& C = composition.Y(speciesi);
    const auto& R = composition.R(speciesi);
    const auto& Deff = composition.Deff(speciesi, U, eps*Saturation);
    const auto& lambda = composition.lambda(speciesi);
    const auto& sourceTerm = composition.sourceTerm(speciesi);

    fvScalarMatrix CEqn
        (
            eps * R * Saturation * fvm::ddt(C)
            + fvm::div(phi, C, "div(phi,C)")
            - fvm::laplacian(eps*Saturation*Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm
            - eps * R * Saturation * fvm::Sp(lambda,C)
        );

    CEqn.solve();

    dCdTmax[speciesi] = max(mag(fvc::ddt(C))).value();

    Info<< "Concentration Min(" << speciesName << ") = " << min(C).value() 
        << " Max(" << speciesName << ") = " << max(C).value() 
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*Saturation*eps).value() 
        << endl;

}
