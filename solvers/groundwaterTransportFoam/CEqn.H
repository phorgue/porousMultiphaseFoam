forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];
    
    auto& C = composition.Y(speciesi);
    const auto& R = composition.R(speciesi);
    const auto& Deff = composition.Deff(speciesi, Utheta, theta);
    const auto& lambda = composition.lambda(speciesi);
    const auto& sourceTerm_tracer = composition.sourceTerm(speciesi);

    fvScalarMatrix CEqn
        (
            R * fvm::ddt(theta,C)
            + fvm::div(phi, C, "div(phi,C)")
            - fvm::laplacian(theta*Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm_tracer
            - R * theta * fvm::Sp(lambda,C)
        );

    CEqn.solve("C");

    dCdTmax[speciesi] = max(mag(fvc::ddt(C))).value();

    Info<< "Concentration Min(" << speciesName << ") = " << min(C).value()
        << " Max(" << speciesName << ") = " << max(C).value()
        << " d" << speciesName << "dTmax = " << dCdTmax[speciesi]
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(C*theta).value()
        << endl;

}
