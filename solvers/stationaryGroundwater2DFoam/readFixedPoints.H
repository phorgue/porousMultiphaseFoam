List<Tuple2<point,scalar> > fixedPotentialList(transportProperties.lookup("fixedPotentialList"));

labelList fixedPotentialIDList(fixedPotentialList.size());
scalarList fixedPotentialValueList(fixedPotentialList.size());
forAll(fixedPotentialList,pointi)
{
    fixedPotentialIDList[pointi] = mesh.findNearestCell(fixedPotentialList[pointi].first());
    fixedPotentialValueList[pointi] = fixedPotentialList[pointi].second();
}

if (MNTpresent)
{
    
}
