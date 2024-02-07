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
    timestepManager& dtm = MDTM.dtManager(speciesi);
    dtm.updateDerivatives();

    Info<< "Concentration: Min(" << speciesName << ") = " << gMin(C.internalField())
        << " Max(" << speciesName << ") = " << gMax(C.internalField())
        << " d" << speciesName << "dTmax = " <<  dtm.dVmax()
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(C*theta).value()
        << " dCmax = " << dtm.dVmax()*runTime.deltaTValue()
        << endl;

}
