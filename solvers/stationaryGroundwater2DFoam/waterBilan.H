{
    //- Compute infiltration removing fixed cells
    Info << "Water bilan (m3/s) : infiltration = " << fvc::domainIntegrate(infiltration).value()*(mesh.nCells()-fixedPotentialIDList.size())/mesh.nCells() << " ; ";

    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi]) << " ; ";
        }
    }

    if (fixedPotentialIDList.size() > 0)
    {
        scalar fixedFlowRate = 0;
        cellFlux = fvc::div(phi*fvc::interpolate(hwater));
        forAll(fixedPotentialIDList,pointi)
        {
            fixedFlowRate += cellFlux.internalField()[fixedPotentialIDList[pointi]]*mesh.V()[fixedPotentialIDList[pointi]];
        }
        Info << " fixed points = " << fixedFlowRate;
    }
    Info << endl;
}
