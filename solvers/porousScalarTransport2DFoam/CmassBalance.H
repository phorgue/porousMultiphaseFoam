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
    Info << " fixed points = " << fvc::domainIntegrate(seepageTerm*C).value()/zScale << endl;

    //- CSV output
    if (CSVactive)
    {
        volTensorField Deff = eps*hwater*DeffModel->Deff();
        if ((!outputEventIsPresent) || outputEvent.currentEventEndTime() == runTime.timeOutputValue())
        {
            massBalanceCSV << runTime.timeName() << " " << fvc::domainIntegrate(R*hwater*C*eps).value()/zScale;
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    scalarField dispersiveFlux = ((Deff.boundaryField()[patchi] & mesh.boundary()[patchi].nf()) & mesh.boundary()[patchi].Sf());
                    scalarField convectiveFlux = phiWater.boundaryField()[patchi]*C.boundaryField()[patchi];
                    massBalanceCSV << " " << sum(dispersiveFlux*fvc::snGrad(C)+convectiveFlux)/zScale;
                }
            }
            massBalanceCSV << " " << fvc::domainIntegrate(seepageTerm*C).value()/zScale << endl;
        }
    }
}
