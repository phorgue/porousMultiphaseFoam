//- Water infiltration event
bool infiltrationEventIsPresent = transportProperties.found("eventFileInfiltration");
word infiltrationEventFileName = transportProperties.lookupOrDefault<word>("eventFileInfiltration","");
infiltrationEventFile infiltrationEvent(infiltrationEventFileName);
if (infiltrationEventIsPresent)
{
    infiltrationEvent.update(runTime.startTime().value());
    infiltrationEvent.storeOldValues();
}

//- Updating infiltration value using event file
if (infiltrationEventIsPresent)
{
    if (infiltrationEvent.isUniform())
    {
        infiltration.primitiveFieldRef() = infiltrationEvent.currentValue(0);
    }
    else
    {
        if (infiltrationEvent.currentValues().size() != mesh.C().size())
        {
            Warning() << "Number of infiltration values does not correspond to the mesh size, first value is used as uniform infiltration" << endl;
            infiltrationEvent.setInfiltrationAsUniform();
            infiltration.primitiveFieldRef() = infiltrationEvent.currentValue(0);
        }
        else
        {
            infiltration.primitiveFieldRef() = infiltrationEvent.currentValues();
        }

    }
}

//- patch event reader for time step managing
forAll(patchEventList,patchEventi)
{
    patchEventList[patchEventi]->update(runTime.startTime().value());
    patchEventList[patchEventi]->storeOldValues();
}

//- Output event reading
bool outputEventIsPresent = runTime.controlDict().found("eventFileOutput");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("eventFileOutput","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent) outputEvent.update(runTime.startTime().value());

scalar eventTimeStep(runTime.controlDict().lookupOrDefault<scalar>("eventTimeStep",0));
if (eventTimeStep > 0)
{
    if (infiltrationEventIsPresent) infiltrationEvent.addIntermediateTimeSteps(eventTimeStep);
}
