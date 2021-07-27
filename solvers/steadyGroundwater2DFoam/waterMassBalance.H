{
    //- Compute infiltration
    Info << "Water mass balance (m3/s) : infiltration = " << (fvc::domainIntegrate(infiltration).value()/zScale) << " ; ";

    seepageTerm = dimensionedScalar("",dimLength/dimTime,0);

    if (fixedPotentialIDList.size() > 0)
    {
        scalar flowOutFixedPoints = 0;
        forAll(fixedPotentialIDList,pointi)
        {
            label currentCell = fixedPotentialIDList[pointi];
            scalar area = mesh.V()[currentCell]/zScale;
            seepageTerm[currentCell] -= cellFlux.internalField()[currentCell];
            flowOutFixedPoints -= cellFlux.internalField()[currentCell]*area;
        }
        Info << " fixed cells = " << flowOutFixedPoints << " ; ";
    }
    if (seepageActive)
    {
        scalar flowOutSeepage = 0;
        forAll(seepageIDList,pointi)
        {
            label currentCell = seepageIDList[pointi];
            scalar area = mesh.V()[currentCell]/zScale;
            seepageTerm[currentCell] -= cellFlux.internalField()[currentCell];
            flowOutSeepage -= cellFlux.internalField()[currentCell]*area;
        }
        Info << " seepage cells = " << flowOutSeepage << " ; ";
    }

    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi])/zScale << " ; ";
        }
    }
    Info << endl;

    //- Residual computation
    scalarField ResiduN = -fvc::laplacian(transmissivity,potential,"laplacian(transmissivity,potential)") + infiltration + seepageTerm;
    forAll(dryCellIDList, celli)  ResiduN[dryCellIDList[celli]] = 0;
    Info << "Potential equation residual = " << gMax(mag(ResiduN)) << endl;
}
