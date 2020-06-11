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
            - fvm::laplacian(eps*hwater*Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm * zScale
            - eps * R * hwater * fvm::Sp(lambda,C)
            - fvm::Sp(seepageTerm,C)
        );

    CEqn.solve("C");

    dCdTmax[speciesi] = max(mag(fvc::ddt(C))).value();
    if (timeScheme == "Euler")
    {
        volScalarField dC2dT2(d2dt2Operator.fvcD2dt2(C));
        dC2dT2max[speciesi] = 0;
        forAll(dC2dT2, celli)
        {
            if(mag(dC2dT2[celli]) > dC2dT2max[speciesi])
            {
                Cmax[speciesi] = mag(C[celli]);
                dC2dT2max[speciesi] = mag(dC2dT2[celli]);
            }
        }
    }
    else
    {
        volScalarField dC3dT3(d3dt3Operator.d3dt3(C));
        dC3dT3max[speciesi] = 0;
        forAll(dC3dT3, celli)
        {
            if(mag(dC3dT3[celli]) > dC3dT3max[speciesi])
            {
                Cmax[speciesi] = mag(C[celli]);
                dC3dT3max[speciesi] = mag(dC3dT3[celli]);
            }
        }
    }
    Info<< "Concentration Min(" << speciesName << ") = " << min(C).value()
        << " Max(" << speciesName << ") = " << max(C).value()
        << " mass(" <<speciesName << ") = " << fvc::domainIntegrate(R*C*hwater*eps).value()/zScale
        << " dCmax = " << dCdTmax[speciesi]*runTime.deltaTValue() << endl;

}
