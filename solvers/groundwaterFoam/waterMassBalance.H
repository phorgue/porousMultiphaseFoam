{
    //- terminal display
    Info << "Water mass balance (m3/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value() << " ; ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;

    //- CSV output
    if (CSVactive)
    {
        if ((!outputEventIsPresent) || outputEvent.currentEventEndTime() == runTime.timeOutputValue())
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

}
