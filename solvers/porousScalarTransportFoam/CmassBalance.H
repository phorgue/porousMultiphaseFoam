forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];
    const auto& C = composition.Y(speciesi);
    const auto& R = composition.R(speciesi);
    const auto sourceTerm = composition.sourceTerm(speciesi);

    //- terminal display
    Info << speciesName << " mass balance (kg/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value() << " ; ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]*C.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;

    //- CSV output
    if (CSVoutput)
    {
        auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];
        
        if ((!outputEventIsPresent) || outputEvent.currentEventEndTime() == runTime.timeOutputValue())
        {
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
