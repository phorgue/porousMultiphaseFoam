{
    //- terminal display
    Info << "C mass balance (kg/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value()/zScale << " ; ";
    forAll(phiWater.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phiWater.boundaryField()[patchi].patch().name() << " = " <<  sum(phiWater.boundaryField()[patchi]*C.boundaryField()[patchi])/zScale << " ; ";
        }
    }
    if (max(seepageTerm).value() > 0)
    {
        Info << " seepage = " << fvc::domainIntegrate(seepageTerm*C).value();
    }
    Info << endl;

    //- CSV output
    if (CSVactive)
    {
        if ((!outputEventIsPresent) || outputEvent.currentEventEndTime() == runTime.timeOutputValue())
        {
            massBalanceCSV << runTime.timeName() << " " << fvc::domainIntegrate(R*hwater*C*eps).value()/zScale;
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    massBalanceCSV << " " << sum(phiWater.boundaryField()[patchi]*C.boundaryField()[patchi])/zScale;
                }
            }
            massBalanceCSV << endl;
        }
    }
}
