//- Fix values if forcingPotential "fixedValue" or "fixedMNT" or "seepage"
if ((forcingPotential.substr(0,5) == "fixed") && (fixedPotentialIDList.size() > 0))
{
    if ( (forcingPotential.substr(0,5) == "fixed") || (forcingPotential == "seepage") )
    {
        potentialEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);
    }
}
else if (forcingPotential == "seepage")
{
    fixedPotentialIDList.clear();
    fixedPotentialValueList.clear();
    forAll(mesh.C(),celli)
    {
        if(potential[celli] >= potentialMNT[celli])
        {
            if (cellFlux[celli] < 0)
            {
                fixedPotentialIDList.append(celli);
                fixedPotentialValueList.append(potentialMNT[celli]);   
            }
        }
    }
    
    Info << "Number of seepage cells = " << fixedPotentialIDList.size() << endl;
    if (fixedPotentialIDList.size() > 0)
    {
        potentialEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);
    }
}
