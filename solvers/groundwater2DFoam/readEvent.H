bool uniformInfiltrationEventIsPresent = transportProperties.found("eventFileInfiltration");
word uniformInfiltrationEventFileName = transportProperties.lookupOrDefault<word>("eventFileInfiltration","");
uniformInfiltrationEventFile uniformInfiltrationEvent(uniformInfiltrationEventFileName);
if (uniformInfiltrationEventIsPresent)
{
    uniformInfiltrationEvent.update(runTime.startTime().value());
    uniformInfiltrationEvent.storeOldValues();
}

bool outputEventIsPresent = runTime.controlDict().found("eventFileOutput");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("eventFileOutput","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent) outputEvent.update(runTime.startTime().value());

scalar eventTimeStep(runTime.controlDict().lookupOrDefault<scalar>("eventTimeStep",0));
if (eventTimeStep > 0)
{
    if (uniformInfiltrationEventIsPresent) uniformInfiltrationEvent.addIntermediateTimeSteps(eventTimeStep);
}

//- Updating infiltration value using event file
if (uniformInfiltrationEventIsPresent)
{
    infiltration.primitiveFieldRef() = uniformInfiltrationEvent.currentValue();
}
