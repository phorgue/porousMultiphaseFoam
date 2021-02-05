//- Fixed potential values
if (fixedPotentialList.size() > 0) potentialEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);

//- Seepage option
if (seepageActive)
{
    seepageIDList.clear();
    seepageValueList.clear();
    forAll(mesh.C(),celli)
    {
        if(potential[celli] >= potentialMNT[celli])
        {
            if (cellFlux[celli] < 0)
            {
                seepageIDList.append(celli);
                seepageValueList.append(potentialMNT[celli]);
            }
        }
    }
    
    Info << "Number of seepage cells = " << seepageIDList.size() << endl;
    if (seepageIDList.size() > 0) potentialEqn.setValues(seepageIDList,seepageValueList);
}
