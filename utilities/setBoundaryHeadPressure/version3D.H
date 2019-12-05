{

    volScalarField h
        (
            IOobject
            (
                "h",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

    //-- Reading patch information
    word patchName = args.opt("patch");
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    fvPatchScalarField& hPatch = h.boundaryFieldRef()[patchID];
    const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
  
    //-- Compute and set up head pressure
    if (args.found("file"))
    {
        //- reading STL informations
        MNTfile fixedPotentialMNT(args.opt("file"));
        Info << nl << "Potential fixed using MNT " << endl;
      
        //- computing local potential
        forAll(hPatch,facei)
        {
            hPatch[facei] = fixedPotentialMNT.interpolate(faces[facei]) - faces[facei].z();
        }
    }  
    else
    {
        scalar potential = args.opt("value",0.);
        Info << nl << "Uniform potential fixed = " << potential << " m " << endl;
        forAll(faces,facei)
        {
            hPatch[facei] = (potential - faces[facei].z());
        }
    }

    h.write();

}
