bool sourceEventIsPresent = transportProperties.found("eventFileTracerSource");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileTracerSource","");
sourceEventFile sourceEvent(sourceEventFileName);
if (sourceEventIsPresent)
{
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

//- determining cell id corresponding to datas
labelList idCoordinates(sourceEvent.ncoordinates());
if (sourceEventIsPresent)
{
    sourceTerm.primitiveFieldRef() = 0;
    forAll(sourceEvent.coordinates(),pointi)
    {
        idCoordinates[pointi] = mesh.findNearestCell(sourceEvent.coordinates()[pointi]);
        sourceTerm[idCoordinates[pointi]] += sourceEvent.currentValue(pointi)/mesh.V()[idCoordinates[pointi]];
    }
}
