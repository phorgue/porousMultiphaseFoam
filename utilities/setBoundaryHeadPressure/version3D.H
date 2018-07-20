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
    word patchName = args.option("patch");
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    fvPatchScalarField& hPatch = h.boundaryFieldRef()[patchID];
    const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
  
    //-- Compute and set up head pressure
    if (args.optionFound("file"))
    {
        //- reading STL informations
        word STLfile = args.option("file");
        triSurfaceMesh potentialSTL(IOobject(STLfile,mesh));
        pointField pPoints = potentialSTL.points();
        Info << nl << "Potential fixed using STL = " << STLfile << endl;
      
        //- computing local potential
        forAll(hPatch,facei)
        {
            scalar xy_distance = GREAT;
            label id_point = -1;
            forAll(pPoints,pointi)
            {
                scalar tmp_dist = Foam::sqrt(pow(pPoints[pointi].x()-faces[facei].x(),2)+pow(pPoints[pointi].y()-faces[facei].y(),2));
                if (tmp_dist < xy_distance)
                {
                    xy_distance = tmp_dist;
                    id_point = pointi;
                }
            }
            hPatch[facei] = pPoints[id_point].z() - faces[facei].z();
        }

    }  
    else
    {
        scalar potential = args.optionLookupOrDefault<scalar>("value",0.);
        Info << nl << "Uniform potential fixed = " << potential << " m " << endl;
        forAll(faces,facei)
        {
            hPatch[facei] = (potential - faces[facei].z());
        }
    }

    h.write();

}
