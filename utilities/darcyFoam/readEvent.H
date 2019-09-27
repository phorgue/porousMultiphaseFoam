bool sourceEventIsPresent = transportProperties.found("eventFileWaterSource");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileWaterSource","");
sourceEventFile sourceEvent(sourceEventFileName);
if (sourceEventIsPresent)
{
    sourceEvent.update(runTime.startTime().value());
    sourceEvent.storeOldValues();
}
