//- Fixed potential values
if (fixedPotentialList.size() > 0) potentialEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);

//- Seepage option
if (seepageActive)
{
    seepageIDList.clear();
    seepageValueList.clear();
    forAll(mesh.C(),celli)
    {
        if(potential[celli] >= potentialDEM[celli])
        {
            if (cellFlux[celli] < 0)
            {
                seepageIDList.append(celli);
                seepageValueList.append(potentialDEM[celli]);
            }
        }
    }

    // Display number of seepage only for serial run (wrong value for parallel
    if (Pstream::nProcs() == 1) Info << "Number of seepage cells = " << seepageIDList.size() << endl;

    if (seepageIDList.size() > 0) potentialEqn.setValues(seepageIDList,seepageValueList);
}
