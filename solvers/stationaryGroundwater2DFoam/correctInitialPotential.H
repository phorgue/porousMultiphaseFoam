{
    if (min(hwater).value() <= 0)
    {
        //- if volume initialization negative
        if (min(hwater.internalField()).value() <= 0)
        {
            Info << "*** WARNING *** : hwater negative in domain (initialized to 0.01) " << endl;

            forAll(mesh.C(),celli)
            {
                hwater[celli] = max(hwater[celli],0.01);
                potential[celli] =  max(
                    z0[celli]+0.01,
                    potential[celli]
                );
            }
        }
        //- if fixed boundary negative
        forAll(mesh.boundary(),patchi)
        {
            if (min(hwater.boundaryField()[patchi]) <= 0)
            {
                if (isA< fixedValueFvPatchField<scalar> >(hwater.boundaryField()[patchi]))
                {
                    Info << "*** WARNING *** : hwater negative fixed (initialized to 0.01) in patch : " << mesh.boundary()[patchi].name() << endl;                 
                    forAll(hwater.boundaryFieldRef()[patchi],facei)
                    {
                        hwater.boundaryFieldRef()[patchi][facei] =  max(0.01,hwater.boundaryFieldRef()[patchi][facei]);
                        potential.boundaryFieldRef()[patchi][facei] =  max(
                            z0.boundaryFieldRef()[patchi][facei]+0.01,
                            potential.boundaryFieldRef()[patchi][facei]
                        );
                    }
                }
            }
        }
        potential.write();
    }
}
