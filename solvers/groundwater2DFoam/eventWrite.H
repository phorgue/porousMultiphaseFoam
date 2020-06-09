//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() <= runTime.timeOutputValue())
    {
        //- compute and write time interpolated fields
        volScalarField potential_interpolated = outputEvent.timeInterpolate(potential, runTime);
        volScalarField hwater_interpolated = outputEvent.timeInterpolate(hwater, runTime);
        volVectorField U_interpolated = outputEvent.timeInterpolate(U, runTime);
        surfaceScalarField phi_interpolated = outputEvent.timeInterpolate(phi, runTime);
        volScalarField infiltration_interpolated = outputEvent.timeInterpolate(infiltration, runTime);
        volScalarField seepageTerm_interpolated = outputEvent.timeInterpolate(seepageTerm, runTime);

        //- CSV output
        if (CSVoutput)
        {
            scalar flowOutFixedPoints_interpolated = outputEvent.timeInterpolate(flowOutFixedPointsPrev,flowOutFixedPoints,runTime);
            waterMassBalanceCSV << outputEvent.currentEventEndTime() << " " <<  (fvc::domainIntegrate(infiltration_interpolated).value()/zScale)*
                (mesh.nCells() - fixedPotentialIDList.size())/mesh.nCells() << " " << flowOutFixedPoints_interpolated << " ";
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    waterMassBalanceCSV << " " << sum(phi_interpolated.boundaryField()[patchi]*hwater_interpolated.boundaryField()[patchi])/zScale << " ";
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
        waterMassBalanceCSV << runTime.timeName() << " " <<  (fvc::domainIntegrate(infiltration).value()/zScale)*(mesh.nCells()-
        fixedPotentialIDList.size())/mesh.nCells() << " " << flowOutFixedPoints << " ";
        forAll(mesh.boundaryMesh(),patchi)
        {
            if (mesh.boundaryMesh()[patchi].type() == "patch")
            {
                waterMassBalanceCSV << " " << sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi])/zScale << " ";
            }
        }
        waterMassBalanceCSV << endl;
    }
}
