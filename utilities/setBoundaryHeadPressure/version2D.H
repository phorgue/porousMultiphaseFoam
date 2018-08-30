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
  
    //-- Compute and set up potential
    
    //- reading STL informations
    word STLfile = args.option("file");
    triSurfaceMesh potentialSTL(IOobject(STLfile,mesh));
    pointField pPoints = potentialSTL.points();
      
    //- computing local potential
    forAll(potentialPatch,facei)
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
        potentialPatch[facei] = pPoints[id_point].z();
    }


    potential.write();
    
}
