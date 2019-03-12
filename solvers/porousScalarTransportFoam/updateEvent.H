if (sourceEventIsPresent)  sourceEvent.update(runTime.timeOutputValue());
if (patchEventIsPresent) patchEvent.update(runTime.timeOutputValue());
if (outputEventIsPresent) outputEvent.update(runTime.timeOutputValue());

if (sourceEventIsPresent)
{
    sourceTerm.primitiveFieldRef() = 0;
    if ((sourceEvent.currentEventStartTime() == runTime.timeOutputValue()) && (isBackwardScheme))
    {
        scalar coefft0_00 = runTime.deltaTValue()/(runTime.deltaTValue() + runTime.deltaT0Value());
        scalar coefftn_0 = 1 + coefft0_00;

        forAll(sourceEvent.coordinates(),pointi)
        {
            sourceTerm[idCoordinates[pointi]] = (coefftn_0*sourceEvent.currentValue(pointi) - coefft0_00*sourceEvent.previousValue(pointi))/mesh.V()[idCoordinates[pointi]];
        }
    }
    else
    {
        forAll(sourceEvent.coordinates(),pointi)
        {
            sourceTerm[idCoordinates[pointi]] = sourceEvent.currentValue(pointi)/mesh.V()[idCoordinates[pointi]];
        }
    }
}
