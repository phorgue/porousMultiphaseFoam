{
    if (min(hwater).value() <= 0)
    {
        //- if volume initialization negative
        if (min(hwater.internalField()).value() <= hwaterMin.value())
        {
            Info << "*** WARNING *** : hwater negative in domain (initialized to hwaterMin) " << endl;

            forAll(mesh.C(),celli)
            {
                hwater[celli] = max(hwater[celli],0.01);
                potential[celli] =  max(
                    z0[celli]+hwaterMin.value(),
                    potential[celli]
                );
            }
        }
        //- if fixed boundary negative
        forAll(mesh.boundary(),patchi)
        {
            if (min(hwater.boundaryField()[patchi]) <= hwaterMin.value())
            {
                if (isA< fixedValueFvPatchField<scalar> >(hwater.boundaryField()[patchi]))
                {
                    Info << "*** WARNING *** : hwater negative fixed (initialized to hwaterMin) in patch : " << mesh.boundary()[patchi].name() << endl;
                    forAll(hwater.boundaryFieldRef()[patchi],facei)
                    {
                        hwater.boundaryFieldRef()[patchi][facei] =  max(0.01,hwater.boundaryFieldRef()[patchi][facei]);
                        potential.boundaryFieldRef()[patchi][facei] =  max(
                            z0.boundaryFieldRef()[patchi][facei]+hwaterMin.value(),
                            potential.boundaryFieldRef()[patchi][facei]
                        );
                    }
                }
            }
        }
        potential.write();
    }
}
