if (infiltrationEventIsPresent)
{
    infiltrationEvent.updateValue(runTime);
    if (infiltrationEvent.isUniform())
    {
        infiltration.primitiveFieldRef() = infiltrationEvent.dtValue(0);
    }
    else
    {
        infiltration.primitiveFieldRef() = infiltrationEvent.dtValues();
    }
}
