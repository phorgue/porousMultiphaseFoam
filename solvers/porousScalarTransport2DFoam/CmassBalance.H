forAll(composition.Y(), speciesi)
{
    const auto& speciesName = composition.species()[speciesi];
    const auto& C = composition.Y(speciesi);
    const auto& sourceTerm = composition.sourceTerm(speciesi);

    //- terminal display
    Info << speciesName << " mass balance (kg/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value() << " ; ";
    forAll(phihwater.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phihwater.boundaryField()[patchi].patch().name() << " = " <<  sum(phihwater.boundaryField()[patchi]*C.boundaryField()[patchi])/zScale << " ; ";
        }
    }
    Info << " fixed points = " << fvc::domainIntegrate(seepageTerm*C).value()/zScale << endl;
}
