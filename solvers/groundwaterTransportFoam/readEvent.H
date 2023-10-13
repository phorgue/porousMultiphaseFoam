//- Water source event reading
word eventFileName_water = transportProperties.lookupOrDefault<word>("waterSourceEventFile","");
bool eventIsPresent_water = transportProperties.found("waterSourceEventFile");
sourceEventFile event_water(eventFileName_water);

if (eventIsPresent_water)
{
    event_water.setTimeScheme(h.name(), mesh);
    event_water.updateIndex(runTime.startTime().value());
    event_water.updateValue(runTime);
    event_water.setFieldDimensions(sourceTerm.dimensions());
}

//- patch event reader for time step managing
forAll(patchEventList,patchEventi)
{
    patchEventList[patchEventi]->updateIndex(runTime.startTime().value());
    patchEventList[patchEventi]->updateValue(runTime);
    patchEventList[patchEventi]->storeOldValues();
}

//- Output event reading
bool outputEventIsPresent = runTime.controlDict().found("outputEventFile");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("outputEventFile","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent)
{
    outputEvent.checkControlDict(runTime);
    outputEvent.updateIndex(runTime.startTime().value());
}
