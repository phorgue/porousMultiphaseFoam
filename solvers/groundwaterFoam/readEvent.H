bool sourceEventIsPresent = transportProperties.found("eventFileSource");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileSource","");
sourceEventFile sourceEvent(sourceEventFileName);
if (sourceEventIsPresent)
{
    sourceEvent.update(runTime.startTime().value());
    sourceEvent.storeOldValues();
}

bool patchEventIsPresent =  transportProperties.found("eventFilePatch");
word patchEventFileName = transportProperties.lookupOrDefault<word>("eventFilePatch","");
patchEventFile patchEvent(patchEventFileName,false);
if (patchEventIsPresent)
{
    patchEvent.update(runTime.startTime().value());
    patchEvent.storeOldValues();
}

bool outputEventIsPresent = runTime.controlDict().found("eventFileOutput");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("eventFileOutput","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent) outputEvent.update(runTime.startTime().value());

scalar eventTimeStep(runTime.controlDict().lookupOrDefault<scalar>("eventTimeStep",0));
if (eventTimeStep > 0)
{
    if (sourceEventIsPresent) sourceEvent.addIntermediateTimeSteps(eventTimeStep);
    if (patchEventIsPresent) patchEvent.addIntermediateTimeSteps(eventTimeStep);
}

//- determining cell id corresponding to datas
labelList idCoordinates(sourceEvent.ncoordinates());
if (sourceEventIsPresent)
{
    sourceTerm.primitiveFieldRef() = 0;
    forAll(sourceEvent.coordinates(),pointi)
    {
        idCoordinates[pointi] = mesh.findNearestCell(sourceEvent.coordinates()[pointi]);
        sourceTerm[idCoordinates[pointi]] = sourceEvent.currentValue(pointi)/mesh.V()[idCoordinates[pointi]];
    }
}
