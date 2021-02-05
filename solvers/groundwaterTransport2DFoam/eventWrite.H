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
        volScalarField waterSourceTerm_interpolated = outputEvent.timeInterpolate(waterSourceTerm, runTime);
        
        //- water CSV output
        if (CSVoutput)
        {
            waterMassBalanceCSV << outputEvent.currentEventEndTime() << " " <<  (fvc::domainIntegrate(infiltration_interpolated).value()/zScale)*
                (mesh.nCells() - fixedPotentialIDList.size())/mesh.nCells() << " ";

            if (waterSourceEventIsPresent) waterMassBalanceCSV << fvc::domainIntegrate(waterSourceTerm_interpolated).value() << " ";
            if (fixedPotentialIDList.size() > 0)
            {
                scalar flowOutFixedPoints_interpolated = outputEvent.timeInterpolate(flowOutFixedPointsPrev,flowOutFixedPoints,runTime);
                waterMassBalanceCSV << flowOutFixedPoints_interpolated << " ";
            }
            if (seepageActive)
            {
                scalar flowOutSeepage_interpolated = outputEvent.timeInterpolate(flowOutSeepagePrev,flowOutSeepage,runTime);
                waterMassBalanceCSV << flowOutSeepage_interpolated << " ";
            }

            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    waterMassBalanceCSV << " " << sum(phi_interpolated.boundaryField()[patchi]*hwater_interpolated.boundaryField()[patchi])/zScale << " ";
                }
            }
            waterMassBalanceCSV << endl;
        }

        //- tracer outputs
        forAll(composition.Y(), speciesi)
        {
            //- compute and write time interpolated fields
            const auto& C = composition.Y(speciesi);
            const auto& R = composition.R(speciesi);
            volScalarField C_interpolated = outputEvent.timeInterpolate(C, runTime);
            volScalarField R_interpolated = outputEvent.timeInterpolate(R, runTime, false);

            //- CSV output at given times only
            if (CSVoutput)
            {
                auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];

                volTensorField Deff(eps*hwater_interpolated*composition.Deff(speciesi));

                CmassBalanceCSV << outputEvent.currentEventEndTime() << " " << fvc::domainIntegrate(R_interpolated*hwater_interpolated*C_interpolated*eps).value()/zScale;
                forAll(mesh.boundaryMesh(),patchi)
                {
                    if (mesh.boundaryMesh()[patchi].type() == "patch")
                    {
                        scalarField dispersiveFlux(((Deff.boundaryField()[patchi] & mesh.boundary()[patchi].nf()) & mesh.boundary()[patchi].Sf()));
                        scalarField convectiveFlux(phi_interpolated.boundaryField()[patchi]*C_interpolated.boundaryField()[patchi]);
                        CmassBalanceCSV << " " << sum(dispersiveFlux*fvc::snGrad(C_interpolated)+convectiveFlux)/zScale;
                    }
                }
                CmassBalanceCSV << " " << fvc::domainIntegrate(seepageTerm*C_interpolated).value()/zScale << endl;
            }
        }

        outputEvent.updateIndex(runTime.timeOutputValue());
    }
}
//- write using openfoam usual runTime
else
{
    //- write C fields at all times
    runTime.write();

    //- write CSV file at all times
    if (CSVoutput)
    {
        //
        waterMassBalanceCSV << runTime.timeName() << " " <<  (fvc::domainIntegrate(infiltration).value()/zScale)*(mesh.nCells()-
        fixedPotentialIDList.size())/mesh.nCells() << " " ;
        if (waterSourceEventIsPresent) waterMassBalanceCSV << fvc::domainIntegrate(waterSourceTerm).value() << " ";
        if (fixedPotentialIDList.size() > 0) Info << flowOutFixedPoints << " ";
        if (seepageActive) Info << flowOutSeepage << " ";

        forAll(mesh.boundaryMesh(),patchi)
        {
            if (mesh.boundaryMesh()[patchi].type() == "patch")
            {
                waterMassBalanceCSV << " " << sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi])/zScale << " ";
            }
        }
        waterMassBalanceCSV << endl;

        forAll(composition.Y(), speciesi)
        {
            auto& CmassBalanceCSV = CmassBalanceCSVs[speciesi];
            const auto& C = composition.Y(speciesi);
            const auto& R = composition.R(speciesi);

            volTensorField Deff(eps*hwater*composition.Deff(speciesi));
            CmassBalanceCSV << runTime.timeName() << " " << fvc::domainIntegrate(R*hwater*C*eps).value()/zScale;
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    scalarField dispersiveFlux(((Deff.boundaryField()[patchi] & mesh.boundary()[patchi].nf()) & mesh.boundary()[patchi].Sf()));
                    scalarField convectiveFlux(phihwater.boundaryField()[patchi]*C.boundaryField()[patchi]);
                    CmassBalanceCSV << " " << sum(dispersiveFlux*fvc::snGrad(C)+convectiveFlux)/zScale;
                }
            }
            CmassBalanceCSV << " " << fvc::domainIntegrate(seepageTerm*C).value()/zScale << endl;
        }
    }
}
