{
    //- terminal display
    Info << "Water mass balance (m3/s) : infiltration = " << fvc::domainIntegrate(infiltration).value()/zScale << " ; ";

    //- computing seepageTerm fields and flowOutFixedPoints
    seepageTerm = dimensionedScalar("",dimLength/dimTime,0);
    flowOutFixedPointsPrev = flowOutFixedPoints;
    flowOutFixedPoints = 0;

    if (fixedPotentialIDList.size() > 0)
    {
        forAll(fixedPotentialIDList,pointi)
        {
            label currentCell = fixedPotentialIDList[pointi];
            scalar area = mesh.V()[currentCell]/zScale;
            seepageTerm[currentCell] -= cellFlux.internalField()[currentCell];
            flowOutFixedPoints -= cellFlux.internalField()[currentCell]*area;
        }
        Info << "fixed points = " << flowOutFixedPoints << " ; ";
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
