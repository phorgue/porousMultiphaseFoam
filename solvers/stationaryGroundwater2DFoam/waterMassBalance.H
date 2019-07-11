{
    //- Compute infiltration removing fixed cells
    Info << "Water mass balance (m3/s) : infiltration = " << (fvc::domainIntegrate(infiltration).value()/zScale)*(mesh.nCells()-fixedPotentialIDList.size())/mesh.nCells() << " ; ";

    seepageTerm = dimensionedScalar("",dimLength/dimTime,0);

    if (fixedPotentialIDList.size() > 0)
    {
        cellFlux = fvc::div(phi*fvc::interpolate(hwater));
        scalar flowOutFixedPoints = 0;
        forAll(fixedPotentialIDList,pointi)
        {
            label currentCell = fixedPotentialIDList[pointi];
            scalar area = mesh.V()[currentCell]/zScale;
            seepageTerm[currentCell] -= cellFlux.internalField()[currentCell]/area;
            flowOutFixedPoints -= cellFlux.internalField()[currentCell]*area;
        }
        Info << " fixed points = " << flowOutFixedPoints;
    }

    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi])/zScale << " ; ";
        }
    }
    Info << endl;
}
