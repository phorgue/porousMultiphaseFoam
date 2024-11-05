//- Fixed potential values
const labelUList fp_ids(fixedPotentialIDList);
const scalarList fp_values(fixedPotentialValueList);
if (fixedPotentialList.size() > 0) potentialEqn.setValues<List>(fp_ids, fp_values);

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

    // Display number of seepage cells
    label nSeepageCells = seepageIDList.size();
    reduce(nSeepageCells, sumOp<label>());
    Info << "Number of seepage cells = " << nSeepageCells << endl;

    if (seepageIDList.size() > 0) potentialEqn.setValues(seepageIDList,seepageValueList);
}
