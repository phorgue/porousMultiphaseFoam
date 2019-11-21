bool sourceEventIsPresent = transportProperties.found("eventFileWaterSource");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileWaterSource","");
sourceEventFile sourceEvent(sourceEventFileName);
if (sourceEventIsPresent)
{
    sourceEvent.setFieldProperties(mesh, sourceTerm.dimensions());

    sourceEvent.update(runTime.startTime().value());
    sourceEvent.storeOldValues();
}

forAll(patchEventList,patchEventi)
{
    patchEventList[patchEventi]->update(runTime.startTime().value());
    patchEventList[patchEventi]->storeOldValues();
}

bool outputEventIsPresent = runTime.controlDict().found("eventFileOutput");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("eventFileOutput","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent) outputEvent.update(runTime.startTime().value());
