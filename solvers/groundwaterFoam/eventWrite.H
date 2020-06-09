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

        //- CSV output
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
    }
}
