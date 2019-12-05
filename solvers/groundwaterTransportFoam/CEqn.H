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

    dCdTmax[speciesi] = max(mag(fvc::ddt(C))).value();
    if (timeScheme == "Euler")
    {
        volScalarField dC2dT2 = d2dt2Operator.fvcD2dt2(C);
        dC2dT2max[speciesi] = 0;
        forAll(dC2dT2, celli)
        {
            if(mag(dC2dT2[celli]) > dC2dT2max[speciesi])
            {
                Cmax[speciesi] = C[celli];
                dC2dT2max[speciesi] = mag(dC2dT2[celli]);
            }
        }
    }
    else
    {
        volScalarField dC3dT3 = d3dt3Operator.d3dt3(C);
        dC3dT3max[speciesi] = 0;
        forAll(dC3dT3, celli)
        {
            if(mag(dC3dT3[celli]) > dC3dT3max[speciesi])
            {
                Cmax[speciesi] = C[celli];
                dC3dT3max[speciesi] = mag(dC3dT3[celli]);
            }
        }
    }

    Info<< "Concentration Min(" << speciesName << ") = " << min(C).value()
        << " Max(" << speciesName << ") = " << max(C).value()
        << " d" << speciesName << "dTmax = " << dCdTmax[speciesi]
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(C*theta).value()
        << " dCmax = " << dCdTmax[speciesi]*runTime.deltaTValue()
        << endl;

}
