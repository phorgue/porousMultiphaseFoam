if (outputEventIsPresent) outputEvent.update(runTime.timeOutputValue());

if (uniformInfiltrationEventIsPresent)
{
    uniformInfiltrationEvent.update(runTime.timeOutputValue());
    infiltration.primitiveFieldRef() = 0;
    if (isBackwardScheme)
    {
        scalar coefft0_00 = runTime.deltaTValue()/(runTime.deltaTValue() + runTime.deltaT0Value());
        scalar coefftn_0 = 1 + coefft0_00;
        infiltration.primitiveFieldRef() = (coefftn_0*uniformInfiltrationEvent.currentValue() - coefft0_00*uniformInfiltrationEvent.oldValue());
    }
    else
    {
        infiltration.primitiveFieldRef() = uniformInfiltrationEvent.currentValue();
    }
}
