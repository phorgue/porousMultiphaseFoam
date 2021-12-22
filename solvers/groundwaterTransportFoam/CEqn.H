//- update dispersion tensor coefficients and source terms
forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
forAll(tracerSourceEventList,tracerSourceEventi) tracerSourceEventList[tracerSourceEventi]->updateValue(runTime);
composition.correct(Utheta, theta);

forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];

    auto& C = composition.Y(speciesi);
    const auto& R = composition.R(speciesi);
    const auto& Deff = composition.Deff(speciesi);
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

    CEqn.solve(mesh.solver("C"));

    dtManagerC[speciesi].updateDerivatives();

    Info<< "Concentration: Min(" << speciesName << ") = " << gMin(C.internalField())
        << " Max(" << speciesName << ") = " << gMax(C.internalField())
        << " d" << speciesName << "dTmax = " << dtManagerC[speciesi].dVmax()
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(C*theta).value()
        << " dCmax = " << dtManagerC[speciesi].dVmax()*runTime.deltaTValue()
        << endl;

}
