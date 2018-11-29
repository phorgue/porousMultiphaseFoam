List<Tuple2<point,scalar> > fixedPotentialList(transportProperties.lookup("fixedPotentialList"));

labelList fixedPotentialIDList(fixedPotentialList.size());
scalarList fixedPotentialValueList(fixedPotentialList.size());
forAll(fixedPotentialList,pointi)
{
    fixedPotentialIDList[pointi] = mesh.findNearestCell(fixedPotentialList[pointi].first());
    fixedPotentialValueList[pointi] = fixedPotentialList[pointi].second();
}

word MNTfileName(transportProperties.lookupOrDefault<word>("MNTfile",""));
if (MNTfileName.size() > 0)
{
    MNTfile fixedPotentialMNT(MNTfileName);
    forAll(fixedPotentialList,pointi)
    {
        fixedPotentialValueList[pointi] = fixedPotentialMNT.interpolate(fixedPotentialList[pointi].first());
    }
}
else
{
    Info <<  " no MNT file"  << endl;
}

if (fixedPotentialList.size() > 0)
{
    Info << " fixed potential positions and values : " << endl;
    forAll(fixedPotentialList,pointi)
    {
        Info << fixedPotentialList[pointi].first() << " " << fixedPotentialValueList[pointi] << endl;
    }
}
