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
  
    //-- Compute and set up head pressure
    if (args.optionFound("MNT"))
    {
        //- reading MNT informations
        MNTfile fixedPotentialMNT(args.option("MNT"));
        Info << nl << "Potential fixed using MNT " << endl;

        //- computing local potential
        forAll(potentialPatch,facei)
        {
            potentialPatch[facei] = fixedPotentialMNT.interpolate(faces[facei]);
        }
    }
    else if (args.optionFound("STL"))
    {
        FatalErrorIn("version2D.H") << "-STL option cannon be used for 2D potential initialization" << abort(FatalError);
    }
    else
    {
        //- uniform potential
        scalar uniformPotential = args.optionLookupOrDefault<scalar>("value",0.);
        Info << nl << "Uniform potential fixed = " << uniformPotential << " m " << endl;
        //- computing and writing local potential
        forAll(potentialPatch,facei)
        {
            potentialPatch[facei] = uniformPotential;
        }
    }

    potential.write();
}
