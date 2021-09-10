//- Read patch name for DEM seepage function                                                                                                                                                              
word patchDEM=transportProperties.lookupOrDefault<word>("patchDEM","none");
label patchDEMID = mesh.boundaryMesh().findPatchID(patchDEM);
labelList topCellID(0);
scalarList distanceToDEM(0);
labelList fixedPotentialIDList(0);
scalarList fixedPotentialValueList(0);

if (patchDEM == "none")
{
    Info << nl << nl << "no DEM patch (no seepage condition)" << endl;
}
else
{
    if (patchDEMID == -1)
    {
        FatalErrorIn("readForcing.H") << "patch for seepage : " << patchDEM << " not found" << abort(FatalError);
    }
    else
    {
        Info << nl << nl << "DEM patch used for seepage = " << patchDEM << " (id=" << patchDEMID << ")" << endl;
        topCellID.resize(mesh.boundaryMesh()[patchDEMID].size());
        topCellID = mesh.boundaryMesh()[patchDEMID].faceCells();
        distanceToDEM.resize(mesh.boundaryMesh()[patchDEMID].size());
        distanceToDEM = mag(mesh.boundary()[patchDEMID].delta()().component(2));
    }
}

