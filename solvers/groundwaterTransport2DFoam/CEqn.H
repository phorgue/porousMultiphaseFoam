//- update dispersion tensor coefficients and source terms
composition.correct(U, eps);

forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];
    
    auto& C = composition.Y(speciesi);
    const auto& R = composition.R(speciesi);
    const auto& Deff = composition.Deff(speciesi);
    const auto& lambda = composition.lambda(speciesi);
    const auto& sourceTerm = composition.sourceTerm(speciesi);

    //- update water flux
    phiWater = phi * fvc::interpolate(hwater);

    fvScalarMatrix CEqn
        (
            eps * R * fvm::ddt(hwater,C)
            + fvm::div(phiWater, C, "div(phi,C)")
            - fvm::laplacian(eps*hwater*Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm * zScale
            - eps * R * hwater * fvm::Sp(lambda,C)
            - fvm::Sp(seepageTerm,C)
        );

    CEqn.solve("C");

    dCdTmax[speciesi] = max(mag(fvc::ddt(C))).value();

    Info<< "Concentration Min(" << speciesName << ") = " << min(C).value()
        << " Max(" << speciesName << ") = " << max(C).value()
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*hwater*eps).value()/zScale << endl;

}
