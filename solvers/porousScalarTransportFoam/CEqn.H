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
            - fvm::laplacian(eps*Saturation*Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm
            - eps * R * Saturation * fvm::Sp(lambda,C)
        );

    CEqn.solve("C");

    dCdTmax[speciesi] = gMax(mag(fvc::ddt(C))());
    if (timeScheme == "Euler")
    {
        volScalarField dC2dT2(d2dt2Operator.fvcD2dt2(C));
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
        volScalarField dC3dT3(d3dt3Operator.d3dt3(C));
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
    Info<< "Concentration: Min(" << speciesName << ") = " << gMin(C.internalField()) 
        << " Max(" << speciesName << ") = " << gMax(C.internalField())
        << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*Saturation*eps).value()
        << " dCmax = " << dCdTmax[speciesi]*runTime.deltaTValue()
        << endl;

}
