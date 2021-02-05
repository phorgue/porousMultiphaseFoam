if (patchMNTID > -1)
{
    fixedPotentialIDList.clear();
    fixedPotentialValueList.clear();
    volScalarField cellFlux(fvc::div(phi));
    forAll(topCellID,celli)
    {
        label currentCell = topCellID[celli];

        if(h[currentCell] >= 0)
        {
            if (cellFlux[currentCell] < 0)
            {
                fixedPotentialIDList.append(currentCell);
                fixedPotentialValueList.append(0);
            }
        }
    }

    Info << "Number of seepage cells = " << fixedPotentialIDList.size() << endl;
    if (fixedPotentialIDList.size() > 0)
    {
        hEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);
    }
}
