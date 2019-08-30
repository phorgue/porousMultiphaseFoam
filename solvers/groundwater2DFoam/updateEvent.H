if (outputEventIsPresent) outputEvent.update(runTime.timeOutputValue());

if (infiltrationEventIsPresent)
{
    infiltrationEvent.update(runTime.timeOutputValue());
    infiltration.primitiveFieldRef() = 0;
    if (isBackwardScheme)
    {
        scalar coefft0_00 = runTime.deltaTValue()/(runTime.deltaTValue() + runTime.deltaT0Value());
        scalar coefftn_0 = 1 + coefft0_00;
        if (infiltrationEvent.isUniform())
        {
            infiltration.primitiveFieldRef() = (coefftn_0*infiltrationEvent.currentValue(0) - coefft0_00*infiltrationEvent.oldValue(0));
        }
        else
        {
            infiltration.primitiveFieldRef() = (coefftn_0*infiltrationEvent.currentValues() - coefft0_00*infiltrationEvent.oldValues());
        }
    }
    else
    {
        if (infiltrationEvent.isUniform())
        {
            infiltration.primitiveFieldRef() = infiltrationEvent.currentValue(0);
        }
        else
        {
            infiltration.primitiveFieldRef() = infiltrationEvent.currentValues();
        }
    }
}
