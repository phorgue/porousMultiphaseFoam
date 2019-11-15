{
    //- water mass balance terminal display
    Info << "Water mass balance (m3/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm_water).value() << " ; ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;

    //- water mass balance CSV output
    if (CSVoutput)
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

    //- tracer mass balance terminal display
    forAll(composition.Y(), speciesi)
    {
        const auto& speciesName = composition.species()[speciesi];
        const auto& C = composition.Y(speciesi);
        const auto& R = composition.R(speciesi);
        const auto& sourceTerm_tracer = composition.sourceTerm(speciesi);

        Info << speciesName << " mass balance (kg/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm_tracer).value() << " ; ";
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
}
