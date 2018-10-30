{
    Info << "*** water bilan ***" << endl;
    forAll(phi.boundaryField(),patchi)
    {
        Info << phi.boundaryField()[patchi].patch().name() << " = " <<  -sum(phi.boundaryField()[patchi]*hwater.boundaryField()[patchi]) << " m3" << endl;
    }
}
