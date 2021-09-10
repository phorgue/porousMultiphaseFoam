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
  
    //-- Compute and set up head pressure
    if (args.found("DEM"))
    {
        //- reading DEM informations
        DEMfile fixedPotentialDEM(args.opt("DEM"));
	scalar offseth = args.optionLookupOrDefault<scalar>("offset",0.);
	if (offseth==0)
	{
	    Info << nl << "Potential fixed using DEM " << endl;
	}
	else
	{
	    Info << nl << "Potential fixed using DEM offseted by " << offseth << endl;
	}

        //- computing local potential
        const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
        forAll(hPatch,facei)
        {
            hPatch[facei] = fixedPotentialDEM.interpolate(faces[facei]) - faces[facei].z() + offseth;
        }
    }
    else if (args.optionFound("STL"))
    {
        //- reading STL informations
        word STLfile = args.opt("STL");
        triSurfaceMesh potentialSTL(IOobject(STLfile,mesh));
        pointField pPoints = potentialSTL.points();
	scalar thresholdh = args.optionLookupOrDefault<scalar>("threshold",0.);
	scalar offseth = args.optionLookupOrDefault<scalar>("offset",0.);
	Info << nl << "Potential fixed using points over " << thresholdh <<  " in STL = " << STLfile << " offseted by " << offseth << endl;
	
        //- computing local potential
        const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
        forAll(hPatch,facei)
        {
            scalar xy_distance = GREAT;
            label id_point = -1;
            forAll(pPoints,pointi)
            {
                if (pPoints[pointi].z()>thresholdh)
		{
                    scalar tmp_dist = Foam::sqrt(pow(pPoints[pointi].x()-faces[facei].x(),2)+pow(pPoints[pointi].y()-faces[facei].y(),2));
                    if (tmp_dist < xy_distance)
		    {
                        xy_distance = tmp_dist;
                        id_point = pointi;
		    }
		}
            }
            hPatch[facei] = pPoints[id_point].z() - faces[facei].z() + offseth;    
        }
    }  
    else
    {
        //- uniform potential
        const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
        scalar potential = args.opt("value",0.);
        Info << nl << "Uniform potential fixed = " << potential << " m " << endl;
        forAll(hPatch,facei)
        {
            hPatch[facei] = (potential - faces[facei].z());
        }
    }

    h.write();

}
