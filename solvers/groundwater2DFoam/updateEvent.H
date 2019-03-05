if (uniformInfiltrationEventIsPresent) uniformInfiltrationEvent.update(runTime.timeOutputValue());
if (outputEventIsPresent) outputEvent.update(runTime.timeOutputValue());

if (uniformInfiltrationEventIsPresent)
{
    if (uniformInfiltrationEvent.currentEventStartTime() == runTime.timeOutputValue())
    {
        infiltration.primitiveFieldRef() = uniformInfiltrationEvent.currentValue();
    }
}
