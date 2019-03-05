if (sourceEventIsPresent) sourceEvent.update(runTime.timeOutputValue());
if (outputEventIsPresent) outputEvent.update(runTime.timeOutputValue());

if (sourceEventIsPresent)
{
    if (sourceEvent.currentEventStartTime() == runTime.timeOutputValue())
    {
        sourceTerm.primitiveFieldRef() = 0;
        forAll(sourceEvent.coordinates(),pointi)
        {
            sourceTerm[idCoordinates[pointi]] = sourceEvent.currentValue(pointi)/mesh.V()[idCoordinates[pointi]];
        }
    }
}
