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

    Info << "Number of seepage cells = " << fixedPotentialIDList.size() << endl;
    if (fixedPotentialIDList.size() > 0)
    {
        hEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);
    }
}
