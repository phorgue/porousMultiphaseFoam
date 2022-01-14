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

    //- update water flux
    phihwater = phi * fvc::interpolate(hwater);
    forAll(dryCellIDList, celli) C[dryCellIDList[celli]] = 0;

    fvScalarMatrix CEqn
        (
            eps * R * fvm::ddt(hwater,C)
            + fvm::div(phihwater, C, "div(phi,C)")
            - fvm::laplacian(eps * hwater * Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm * zScale
            - eps * R * hwater * fvm::Sp(lambda,C)
            - fvm::Sp(seepageTerm,C)
        );

    CEqn.solve(mesh.solver("C"));

    dtManagerC[speciesi].updateDerivatives();

    Info<< "Concentration Min(" << speciesName << ") = " << min(C).value()
        << " Max(" << speciesName << ") = " << max(C).value()
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*hwater*eps).value()/zScale
        << " dCmax = " << dtManagerC[speciesi].dVmax()*runTime.deltaTValue() << endl;

}
