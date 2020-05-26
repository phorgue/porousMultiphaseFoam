forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];
    const auto& C = composition.Y(speciesi);
    const auto& sourceTerm = composition.sourceTerm(speciesi);

    //- terminal display
    Info << speciesName << " mass balance (kg/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value() << " ; ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  sum(phi.boundaryField()[patchi]*C.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;

}
