//- update dispersion tensor coefficients and source terms
forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
forAll(tracerSourceEventList,tracerSourceEventi) tracerSourceEventList[tracerSourceEventi]->updateValue(runTime);

//- Correct pmTransportModel + dispersion for classical porosity
pmTransportModel->solveTransport(Utheta, phi, theta);

//- display informations
forAll(composition.Y(), speciesi)
{
    const volScalarField& C = composition.Y(speciesi);
    const word& speciesName = composition.species()[speciesi];
    dtManagerC[speciesi].updateDerivatives();

    Info<< "Concentration: Min(" << speciesName << ") = " << gMin(C.internalField())
        << " Max(" << speciesName << ") = " << gMax(C.internalField())
        << " d" << speciesName << "dTmax = " << dtManagerC[speciesi].dVmax()
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(C*theta).value()
        << " dCmax = " << dtManagerC[speciesi].dVmax()*runTime.deltaTValue()
        << endl;

}
