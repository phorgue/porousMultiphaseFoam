//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() < runTime.timeOutputValue())
    {
        scalar interpolateFactor = (runTime.timeOutputValue()-outputEvent.currentEventEndTime())/ runTime.deltaTValue();
        scalar timeOutputBackup = runTime.timeOutputValue();
        runTime.setTime(outputEvent.currentEventEndTime(),runTime.timeIndex());

        forAll(composition.Y(), speciesi)
        {
            //- Write C fields at given times only
            const auto& speciesName = composition.species()[speciesi];
            const auto& C = composition.Y(speciesi);
            const auto& R = composition.R(speciesi);

            volScalarField C_interpolated
                (
                    IOobject
                    (
                        speciesName,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    C
                );
            C_interpolated = interpolateFactor*C+(1.0-interpolateFactor)*C.oldTime();
            C_interpolated.write();

            volScalarField R_interpolated = interpolateFactor*R+(1.0-interpolateFactor)*R.oldTime();

            //- CSV output at given times only
            if (CSVoutput)
            {
                auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];
                Info << "passe ici" << endl;
                CmassBalanceCSV << runTime.timeName() << " " << fvc::domainIntegrate(R_interpolated*C_interpolated*Saturation*eps).value();
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

        runTime.setTime(timeOutputBackup,runTime.timeIndex());
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
