if (patchDEMID > -1)
{
    fixedPotentialIDList.clear();
    fixedPotentialValueList.clear();
    volScalarField cellFlux(fvc::div(phi));
    forAll(topCellID,celli)
    {
        label currentCell = topCellID[celli];

        if(h[currentCell] >= distanceToDEM[celli])
        {
            if (cellFlux[currentCell] < 0)
            {
                fixedPotentialIDList.append(currentCell);
                fixedPotentialValueList.append(distanceToDEM[celli]);
            }
        }
    }

    // Display number of seepage only for serial run (wrong value for parallel
    if (Pstream::nProcs() == 1) Info << "Number of seepage cells = " << fixedPotentialIDList.size() << endl;

    if (fixedPotentialIDList.size() > 0) hEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);
}
