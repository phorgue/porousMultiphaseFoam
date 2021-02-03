{
    //- water mass balance terminal display
    Info << "Water mass balance (m3/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value() << " ; ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;

    //- tracer mass balance terminal display
    forAll(composition.Y(), speciesi)
    {
        const auto& speciesName = composition.species()[speciesi];
        const auto& C = composition.Y(speciesi);
        const auto& sourceTerm_tracer = composition.sourceTerm(speciesi);

        Info << speciesName << " mass balance (kg/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm_tracer).value() << " ; ";
        forAll(phi.boundaryField(),patchi)
        {
            if (mesh.boundaryMesh()[patchi].type() == "patch")
            {
                Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]*C.boundaryField()[patchi]) << " ; ";
            }
        }
        Info << endl;
    }
}
