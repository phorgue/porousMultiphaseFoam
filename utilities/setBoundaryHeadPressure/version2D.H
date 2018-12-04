{

    volScalarField potential
        (
            IOobject
            (
                "potential",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
 
    //-- Reading patch information
    word patchName = args.option("patch");
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    fvPatchScalarField& potentialPatch = potential.boundaryFieldRef()[patchID];
    const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
  
    //- reading MNT informations
    MNTfile fixedPotentialMNT(args.option("file"));

    //- computing and writing local potential
    forAll(potentialPatch,facei)
    {
        potentialPatch[facei] = fixedPotentialMNT.interpolate(faces[facei]);
    }
    potential.write();
    
}
