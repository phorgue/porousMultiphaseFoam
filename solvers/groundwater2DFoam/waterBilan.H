{
    scalar waterBilan = fvc::domainIntegrate(infiltration).value();
    Info << "Water bilan : ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
        scalar patchFlux = sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi]);
        waterBilan -= patchFlux;
        Info << phi.boundaryField()[patchi].patch().name() << " = " <<  -patchFlux << " m3 ; ";
        }
    }
    Info << "TOTAL = " << waterBilan << " m3" << endl;
    if (mag(waterBilan) < stationaryTolerance)
    {
        Info << nl << "Stationary state reached : write and end run" << endl;
        runTime.writeAndEnd();
    }
}
