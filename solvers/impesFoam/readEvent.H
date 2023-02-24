bool sourceEventIsPresent = transportProperties.found("sourceEventFileWater");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("sourceEventFileWater","");
sourceEventFile sourceEvent(sourceEventFileName);
if (sourceEventIsPresent)
{
    sourceEvent.setTimeScheme(Sb.name(), mesh);
    sourceEvent.setFieldDimensions(sourceTerm.dimensions());

    sourceEvent.updateIndex(runTime.startTime().value());
    sourceEvent.updateValue(runTime);
    sourceEvent.storeOldValues();
}

forAll(patchEventList,patchEventi)
{
    patchEventList[patchEventi]->updateIndex(runTime.startTime().value());
    patchEventList[patchEventi]->updateValue(runTime);
    patchEventList[patchEventi]->storeOldValues();
}

bool outputEventIsPresent = runTime.controlDict().found("outputEventFile");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("outputEventFile","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent)
{
    outputEvent.checkControlDict(runTime);
    outputEvent.updateIndex(runTime.startTime().value());
}
