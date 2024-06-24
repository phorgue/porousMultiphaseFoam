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
    if (args.found("DEM"))
    {
        //- reading DEM informations
        DEMfile fixedPotentialDEM(args.opt("DEM"));
	scalar offsetPotential = args.getOrDefault<scalar>("offset",0.);
	if (offsetPotential==0)
	{
	    Info << nl << "Potential fixed using DEM" << endl;
	}
	else
	{
	    Info << nl << "Potential fixed using DEM offseted by " << offsetPotential << endl;
	}

        //- computing local potential
        const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
        forAll(potentialPatch,facei)
        {
            potentialPatch[facei] = fixedPotentialDEM.interpolate(faces[facei]) + offsetPotential;
        }
    }
    else if (args.found("STL"))
    {
        //- reading STL informations
        word STLfile = args.opt("STL");
        triSurfaceMesh potentialSTL(IOobject(STLfile,mesh));
        pointField pPoints( potentialSTL.points() );
	scalar thresholdPotential = args.getOrDefault<scalar>("threshold",0.);
	scalar offsetPotential = args.getOrDefault<scalar>("offset",0.);
        Info << nl << "Potential fixed using points over " << thresholdPotential <<  " in STL = " << STLfile << " offseted by " << offsetPotential << endl;
        const vectorField faces = mesh.boundary()[patchID].patch().faceCentres();
	
        //- computing local potential
        forAll(potentialPatch,facei)
        {
            scalar xy_distance = GREAT;
            label id_point = -1;
            forAll(pPoints,pointi)
            {
                if (pPoints[pointi].z()>thresholdPotential)
		{
                    scalar tmp_dist = Foam::sqrt(pow(pPoints[pointi].x()-faces[facei].x(),2)+pow(pPoints[pointi].y()-faces[facei].y(),2));
                    if (tmp_dist < xy_distance)
		    {
                        xy_distance = tmp_dist;
                        id_point = pointi;
		    }
		}
            }
            potentialPatch[facei] = pPoints[id_point].z() + offsetPotential;    
        }
    }  
    else if (args.found("value"))
    {
        //- uniform potential
        scalar uniformPotential = args.getOrDefault<scalar>("value",0.);
        Info << nl << "Uniform potential fixed = " << uniformPotential << " m " << endl;
        //- computing and writing local potential
        forAll(potentialPatch,facei)
        {
            potentialPatch[facei] = uniformPotential;
        }
    }

    potential.write();
}
