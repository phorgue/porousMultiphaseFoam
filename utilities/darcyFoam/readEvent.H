bool sourceEventIsPresent = transportProperties.found("eventFileWaterSource");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileWaterSource","");
sourceEventFile sourceEvent(sourceEventFileName);
if (sourceEventIsPresent)
{
    sourceEvent.updateIndex(runTime.startTime().value());
    sourceEvent.updateValue(runTime);
    sourceEvent.storeOldValues();
}
