{
    Info << "Water bilan (m3/s) : infiltration = " << fvc::domainIntegrate(infiltration).value() << " ; ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  -sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi]) << " ; ";
        }
    }

    if (fixedPotentialList.size() > 0)
    {
        scalar fixedFlowRate = 0;
        volScalarField contError(fvc::div(phi*fvc::interpolate(hwater)));
        forAll(fixedPotentialIDList,pointi)
        {
            fixedFlowRate += contError.internalField()[fixedPotentialIDList[pointi]]*mesh.V()[fixedPotentialIDList[pointi]];
        }
        Info << " fixed points = " << fixedFlowRate;
    }
    Info << endl;
}
