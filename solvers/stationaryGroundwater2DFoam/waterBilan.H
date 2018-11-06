{
    Info << "Water bilan (m3/s) : ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  -sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;
}
