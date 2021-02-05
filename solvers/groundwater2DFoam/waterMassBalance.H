{
    //- terminal display
    Info << "Water mass balance (m3/s) : infiltration = " << fvc::domainIntegrate(infiltration).value()/zScale << " ; ";

    //- compute source term
    if (waterSourceEventIsPresent) Info << "source = " << fvc::domainIntegrate(waterSourceTerm).value() << " ; ";

    //- computing seepageTerm fields and global flow in/out (for fixed/seepage)
    seepageTerm = dimensionedScalar("",dimLength/dimTime,0);

    flowOutFixedPointsPrev = flowOutFixedPoints;
    flowOutFixedPoints = 0;
    if (fixedPotentialIDList.size() > 0)
    {
        forAll(fixedPotentialIDList,pointi)
        {
            label currentCell = fixedPotentialIDList[pointi];
            scalar area = mesh.V()[currentCell]/zScale;
            flowOutFixedPoints -= cellFlux.internalField()[currentCell]*area;
        }
        Info << "fixed points = " << flowOutFixedPoints << " ; ";
    }

    flowOutSeepagePrev = flowOutSeepage;
    flowOutSeepage = 0;
    if (seepageActive)
    {
        forAll(seepageIDList,pointi)
        {
            label currentCell = seepageIDList[pointi];
            scalar area = mesh.V()[currentCell]/zScale;
            seepageTerm[currentCell] -= cellFlux.internalField()[currentCell];
            flowOutSeepage -= cellFlux.internalField()[currentCell]*area;
        }
        Info << " seepage cells = " << flowOutSeepage << " ; ";
    }

    //- compute/display flow at boundaries
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi])/zScale << " ; ";
        }
    }
    Info << endl;
}
