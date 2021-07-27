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
    word patchName = args.opt("patch");
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    fvPatchScalarField& potentialPatch = potential.boundaryFieldRef()[patchID];

    //-- Compute and set up head pressure
    if (args.optionFound("DEM"))
    {
        //- reading DEM informations
        DEMfile fixedPotentialDEM(args.opt("DEM"));
        Info << nl << "Potential fixed using DEM " << endl;

        //- computing local potential
        const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
        forAll(potentialPatch,facei)
        {
            potentialPatch[facei] = fixedPotentialDEM.interpolate(faces[facei]);
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
