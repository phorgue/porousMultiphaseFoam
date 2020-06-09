//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() <= runTime.timeOutputValue())
    {
        forAll(composition.Y(), speciesi)
        {
            //- compute and write time interpolated fields
            const auto& C = composition.Y(speciesi);
            const auto& R = composition.R(speciesi);
            volScalarField C_interpolated = outputEvent.timeInterpolate(C, runTime);
            volScalarField R_interpolated = outputEvent.timeInterpolate(R, runTime,false);

            //- CSV output at given times only
            if (CSVoutput)
            {
                auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];
                CmassBalanceCSV << outputEvent.currentEventEndTime() << " " << fvc::domainIntegrate(R_interpolated*C_interpolated*Saturation*eps).value();
                forAll(mesh.boundaryMesh(),patchi)
                {
                    if (mesh.boundaryMesh()[patchi].type() == "patch")
                    {
                        CmassBalanceCSV << " " << sum(phi.boundaryField()[patchi]*C_interpolated.boundaryField()[patchi]);
                    }
                }
                CmassBalanceCSV << endl;

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
            const auto& C = composition.Y(speciesi);
            const auto& R = composition.R(speciesi);
            auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];

            CmassBalanceCSV << runTime.timeName() << " " << fvc::domainIntegrate(R*C*Saturation*eps).value();
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    CmassBalanceCSV << " " << sum(phi.boundaryField()[patchi]*C.boundaryField()[patchi]);
                }
            }
            CmassBalanceCSV << endl;
        }
    }
}
