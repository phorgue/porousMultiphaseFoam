//- update dispersion tensor coefficients and source terms
forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
forAll(sourceEventList,sourceEventi) sourceEventList[sourceEventi]->updateValue(runTime);
composition.correct(U, Saturation, eps);

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
            eps * R * Saturation * fvm::ddt(C)
            + fvm::div(phi, C, "div(phi,C)")
            -  fvm::laplacian(eps * Saturation * Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm
            - eps * R * Saturation * fvm::Sp(lambda,C)
        );

    CEqn.solve(mesh.solver("C"));

    dtManager[speciesi].updateDerivatives();

    Info<< "Concentration: Min(" << speciesName << ") = " << gMin(C.internalField()) 
        << " Max(" << speciesName << ") = " << gMax(C.internalField())
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*Saturation*eps).value()
        << " dCmax = " << dtManager[speciesi].dVmax()*runTime.deltaTValue()
        << endl;

}
