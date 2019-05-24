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
            massBalanceCSV << runTime.timeName() << " ";
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    massBalanceCSV << " " << sum(phi.boundaryField()[patchi]);
                }
            }
            massBalanceCSV << endl;
        }
    }

}
