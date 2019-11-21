if (infiltrationEventIsPresent)
{
    if (infiltrationEvent.isUniform())
    {
        infiltration.primitiveFieldRef() = infiltrationEvent.dtValue(0, runTime, mesh.ddtScheme("source"));
    }
    else
    {
        infiltration.primitiveFieldRef() = infiltrationEvent.dtValues(runTime, mesh.ddtScheme("source"));
    }
}
