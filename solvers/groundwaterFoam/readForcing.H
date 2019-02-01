//- Read patch name for MNT seepage function
word patchMNT=transportProperties.lookupOrDefault<word>("patchMNT","none");
label patchMNTID = mesh.boundaryMesh().findPatchID(patchMNT);
labelList topCellID(0);
labelList fixedPotentialIDList(0);
scalarList fixedPotentialValueList(0);

if (patchMNT == "none")
{
    Info << nl << nl << "no MNT patch (no seepage condition)" << endl;
}
else
{
    if (patchMNTID == -1)
    {
        FatalErrorIn("readForcing.H") << "patch for seepage : " << patchMNT << " not found" << abort(FatalError);
    }
    else
    {
        Info << nl << nl << "MNT patch used for seepage = " << patchMNT << " (id=" << patchMNTID << ")" << endl;
        topCellID.resize(mesh.boundaryMesh()[patchMNTID].size());
        topCellID = mesh.boundaryMesh()[patchMNTID].faceCells();
    }
}
