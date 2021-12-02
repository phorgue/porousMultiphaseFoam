//- update dispersion tensor coefficients and source terms
forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
forAll(sourceEventList,sourceEventi) sourceEventList[sourceEventi]->updateValue(runTime);
composition.correct(U, eps);

forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];

    auto& C = composition.Y(speciesi);
    const auto& R = composition.R(speciesi);
    const auto& Deff = composition.Deff(speciesi);
    const auto& lambda = composition.lambda(speciesi);
    const auto& sourceTerm = composition.sourceTerm(speciesi);

    fvScalarMatrix CEqn
        (
            eps * R * hwater * fvm::ddt(C)
            + fvm::div(phihwater, C, "div(phi,C)")
            - fvm::laplacian(eps * hwater * Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm * zScale
            - eps * R * hwater * fvm::Sp(lambda,C)
            - fvm::Sp(seepageTerm,C)
        );

    CEqn.solve(mesh.solver("C"));

    dtManager[speciesi].updateDerivatives();

    Info<< "Concentration: Min(" << speciesName << ") = " << gMin(C.internalField())
        << " Max(" << speciesName << ") = " << gMax(C.internalField())
        << " mass(" <<speciesName << ") = " << fvc::domainIntegrate(R*C*hwater*eps).value()/zScale
        << " dCmax = " << dtManager[speciesi].dVmax()*runTime.deltaTValue() << endl;

}
