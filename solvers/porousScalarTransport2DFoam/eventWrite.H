//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() <= runTime.timeOutputValue())
    {
        volScalarField hwater_interpolated = outputEvent.timeInterpolate(hwater, runTime);
        surfaceScalarField phihwater_interpolated =  outputEvent.timeInterpolate(phihwater, runTime);
        forAll(composition.Y(), speciesi)
        {
            //- compute and write time interpolated fields
            const auto& C = composition.Y(speciesi);
            const auto& R = composition.R(speciesi);
            volScalarField C_interpolated = outputEvent.timeInterpolate(C, runTime);
            volScalarField R_interpolated = outputEvent.timeInterpolate(R, runTime);

            //- CSV output at given times only
            if (CSVoutput)
            {
                auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];

                volTensorField Deff(eps*hwater_interpolated*composition.Deff(speciesi));

                CmassBalanceCSV << outputEvent.currentEventEndTime() << " " << fvc::domainIntegrate(R_interpolated*hwater_interpolated*C_interpolated*eps).value()/zScale;
                forAll(mesh.boundaryMesh(),patchi)
                {
                    if (mesh.boundaryMesh()[patchi].type() == "patch")
                    {
                        scalarField dispersiveFlux(((Deff.boundaryField()[patchi] & mesh.boundary()[patchi].nf()) & mesh.boundary()[patchi].Sf()));
                        scalarField convectiveFlux(phihwater_interpolated.boundaryField()[patchi]*C_interpolated.boundaryField()[patchi]);
                        CmassBalanceCSV << " " << sum(dispersiveFlux*fvc::snGrad(C_interpolated)+convectiveFlux)/zScale;
                    }
                }
                CmassBalanceCSV << " " << fvc::domainIntegrate(seepageTerm*C_interpolated).value()/zScale << endl;
            }
        }

        outputEvent.updateIndex(runTime.timeOutputValue());
    }
}
//- write using openfoam usual runTime
else
{
    //- write C fields at all times
    runTime.write();

    //- write CSV file at all times
    if (CSVoutput)
    {
        forAll(composition.Y(), speciesi)
        {
            auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];
            const auto& C = composition.Y(speciesi);
            const auto& R = composition.R(speciesi);

            volTensorField Deff(eps*hwater*composition.Deff(speciesi));
            CmassBalanceCSV << runTime.timeName() << " " << fvc::domainIntegrate(R*hwater*C*eps).value()/zScale;
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    scalarField dispersiveFlux(((Deff.boundaryField()[patchi] & mesh.boundary()[patchi].nf()) & mesh.boundary()[patchi].Sf()));
                    scalarField convectiveFlux(phihwater.boundaryField()[patchi]*C.boundaryField()[patchi]);
                    CmassBalanceCSV << " " << sum(dispersiveFlux*fvc::snGrad(C)+convectiveFlux)/zScale;
                }
            }
            CmassBalanceCSV << " " << fvc::domainIntegrate(seepageTerm*C).value()/zScale << endl;
        }
    }
}
