List<Tuple2<point,scalar> > fixedPotentialList(transportProperties.lookupOrDefault
("fixedPotentialList",List<Tuple2<point,scalar> >())
);;

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
    Info << nl << "No MNT file"  << endl;
}

if (fixedPotentialList.size() > 0)
{
    Info << nl << "Fixed potential positions and values are " << nl << "{";
    forAll(fixedPotentialList,pointi)
    {
        Info << nl << "  " << fixedPotentialList[pointi].first() << " " << fixedPotentialValueList[pointi];
    }
    Info << nl << "}" << endl;
}
{
    Info << nl << "No fixed potential" << endl;
}
