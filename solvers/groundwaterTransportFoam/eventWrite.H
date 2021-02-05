//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() <= runTime.timeOutputValue())
    {
        //- compute and write time interpolated fields
        volScalarField h_interpolated = outputEvent.timeInterpolate(h, runTime);
        volScalarField theta_interpolated = outputEvent.timeInterpolate(theta, runTime);
        volVectorField Utheta_interpolated = outputEvent.timeInterpolate(Utheta, runTime);
        surfaceScalarField phi_interpolated = outputEvent.timeInterpolate(phi, runTime);

        //- water CSV output
        if (CSVoutput)
        {
            waterMassBalanceCSV << outputEvent.currentEventEndTime() << " ";
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    waterMassBalanceCSV << " " << sum(phi_interpolated.boundaryField()[patchi]);
                }
            }
            waterMassBalanceCSV << endl;
        }

        //- tracer outputs
        forAll(composition.Y(), speciesi)
        {
            //- compute and write time interpolated fields
            const auto& C = composition.Y(speciesi);
            const auto& R = composition.R(speciesi);
            volScalarField C_interpolated = outputEvent.timeInterpolate(C, runTime);
            volScalarField R_interpolated = outputEvent.timeInterpolate(R, runTime, false);

            //- CSV output at given times only
            if (CSVoutput)
            {
                auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];
                CmassBalanceCSV << outputEvent.currentEventEndTime() << " " << fvc::domainIntegrate(R_interpolated*C_interpolated*theta_interpolated).value();
                forAll(mesh.boundaryMesh(),patchi)
                {
                    if (mesh.boundaryMesh()[patchi].type() == "patch")
                    {
                        CmassBalanceCSV << " " << sum(phi_interpolated.boundaryField()[patchi]*C_interpolated.boundaryField()[patchi]);
                    }
                }
                CmassBalanceCSV << endl;
            }
        }

        outputEvent.updateIndex(runTime.timeOutputValue());
    }
}
else
{
    //- write fields using usual openfoam rules
    runTime.write();

    //- write CSV at all times
    if (CSVoutput)
    {
            waterMassBalanceCSV << runTime.timeName() << " ";
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    waterMassBalanceCSV << " " << sum(phi.boundaryField()[patchi]);
                }
            }
            waterMassBalanceCSV << endl;

            forAll(composition.Y(), speciesi)
            {
                const auto& C = composition.Y(speciesi);
                const auto& R = composition.R(speciesi);
                auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];

                CmassBalanceCSV << runTime.timeName() << " " << fvc::domainIntegrate(R*C*theta).value();
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
