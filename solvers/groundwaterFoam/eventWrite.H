//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() < runTime.timeOutputValue())
    {
        scalar interpolateFactor = (runTime.timeOutputValue()-outputEvent.currentEventEndTime())/ runTime.deltaTValue();
        scalar timeOutputBackup = runTime.timeOutputValue();
        runTime.setTime(outputEvent.currentEventEndTime(),runTime.timeIndex());
        volScalarField h_save
            (
                IOobject
                (
                    "h",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                h
            );
        h_save = interpolateFactor*h+(1.0-interpolateFactor)*h.oldTime();
        h_save.write();
        volScalarField theta_save
        
            (
                IOobject
                (
                    "theta",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                theta
            );
        theta_save = interpolateFactor*theta+(1.0-interpolateFactor)*theta.oldTime();
        theta_save.write();
        volVectorField Utheta_save
            (
                IOobject
                (
                    "Utheta",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Utheta
            );
        Utheta_save = interpolateFactor*Utheta+(1.0-interpolateFactor)*Utheta.oldTime();
        Utheta_save.write();
        surfaceScalarField phi_save
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi
            );
        phi_save = interpolateFactor*phi+(1.0-interpolateFactor)*phi.oldTime();
        phi_save.write();

        //- CSV output
        if (CSVoutput)
        {
            waterMassBalanceCSV << runTime.timeName() << " ";
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    waterMassBalanceCSV << " " << sum(phi_save.boundaryField()[patchi]);
                }
            }
            waterMassBalanceCSV << endl;
        }

        runTime.setTime(timeOutputBackup,runTime.timeIndex());
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
