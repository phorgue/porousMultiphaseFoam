if (patchDEMID > -1)
{
    volScalarField cellFlux(fvc::div(phi));
    labelList fixedPotentialIDListPrev(fixedPotentialIDList);
    fixedPotentialIDList.clear();
    fixedPotentialValueList.clear();

    forAll(topCellID,celli)
    {
        label currentCell = topCellID[celli];

        if(h[currentCell] >= distanceToDEM[celli])
        {
            // if(cellFlux[currentCell] > 0)
            // {
            //     bool isFound = false;
            //     forAll(fixedPotentialIDListPrev,pointi)
            //     {
            //         if (fixedPotentialIDListPrev[pointi] == currentCell)
            //         {
            //             isFound = true;
            //         }
            //     }
            //     if (!isFound)
            //     {
            //         fixedPotentialIDList.append(currentCell);
            //         fixedPotentialValueList.append(distanceToDEM[celli]);
            //     }
            // }
            // else
            if (cellFlux[currentCell] < 0)
            {
                fixedPotentialIDList.append(currentCell);
                fixedPotentialValueList.append(distanceToDEM[celli]);
            }
        }
    }

    // Display number of seepage only for serial run (wrong value for parallel
    label nSeepageCells = fixedPotentialIDList.size();
    reduce(nSeepageCells, sumOp<label>());
    Info << "Number of seepage cells = " << nSeepageCells << endl;

    if (fixedPotentialIDList.size() > 0) hEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);
}
