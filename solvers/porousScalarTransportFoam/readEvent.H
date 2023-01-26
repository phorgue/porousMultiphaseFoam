forAll(tracerSourceEventList,sourceEventi)
{
	tracerSourceEventList[sourceEventi]->updateIndex(runTime.startTime().value());
	tracerSourceEventList[sourceEventi]->updateValue(runTime);
	tracerSourceEventList[sourceEventi]->storeOldValues();
}

forAll(patchEventList,patchEventi)
{
    patchEventList[patchEventi]->updateIndex(runTime.startTime().value());
    patchEventList[patchEventi]->updateValue(runTime);
    patchEventList[patchEventi]->storeOldValues();
}

bool outputEventIsPresent = runTime.controlDict().found("eventFileOutput");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("eventFileOutput","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent)
{
    outputEvent.checkControlDict(runTime);
    outputEvent.updateIndex(runTime.startTime().value());
}
