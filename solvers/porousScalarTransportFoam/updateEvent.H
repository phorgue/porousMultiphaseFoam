if (sourceEventIsPresent)
{
    if (runTime.timeOutputValue() >= sourceEvent.dates()[currentSourceEvent+1])
    {
        sourceTerm = dimensionedScalar("",dimMass/(dimVol*dimTime),0);
        currentSourceEvent++;
        forAll(sourceEvent.coordinates(),pointi)
        {
            sourceTerm[idCoordinates[pointi]] = sourceEvent.datas()[currentSourceEvent][pointi];
        }
    }
}

if (patchEventIsPresent)
{
    if (runTime.timeOutputValue() >= patchEvent.dates()[currentPatchEvent+1])
    {
        currentPatchEvent++;
    }
}
