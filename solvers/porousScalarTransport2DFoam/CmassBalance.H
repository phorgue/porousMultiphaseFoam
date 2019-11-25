forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];
    const auto& C = composition.Y(speciesi);
    const auto& R = composition.R(speciesi);
    const auto& sourceTerm = composition.sourceTerm(speciesi);

    //- terminal display
    Info << speciesName << " mass balance (kg/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value()/zScale << " ; ";
    forAll(phiWater.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phiWater.boundaryField()[patchi].patch().name() << " = " <<  sum(phiWater.boundaryField()[patchi]*C.boundaryField()[patchi])/zScale << " ; ";
        }
    }
    Info << " fixed points = " << fvc::domainIntegrate(seepageTerm*C).value()/zScale << endl;

    //- CSV output
    if (CSVoutput)
    {
        auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];
        
        volTensorField Deff = eps*hwater*composition.Deff(speciesi);
        if ((!outputEventIsPresent) || outputEvent.currentEventEndTime() == runTime.timeOutputValue())
        {
            CmassBalanceCSV << runTime.timeName() << " " << fvc::domainIntegrate(R*hwater*C*eps).value()/zScale;
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    scalarField dispersiveFlux = ((Deff.boundaryField()[patchi] & mesh.boundary()[patchi].nf()) & mesh.boundary()[patchi].Sf());
                    scalarField convectiveFlux = phiWater.boundaryField()[patchi]*C.boundaryField()[patchi];
                    CmassBalanceCSV << " " << sum(dispersiveFlux*fvc::snGrad(C)+convectiveFlux)/zScale;
                }
            }
            CmassBalanceCSV << " " << fvc::domainIntegrate(seepageTerm*C).value()/zScale << endl;
        }
    }
}
