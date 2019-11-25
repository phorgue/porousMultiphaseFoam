//- Water source event reading
word eventFileName_water = transportProperties.lookupOrDefault<word>("eventFileWaterSource","");
bool eventIsPresent_water = transportProperties.found("eventFileWaterSource");
sourceEventFile event_water(eventFileName_water);

volScalarField sourceTerm_water
(
    IOobject
    (
        "sourceTerm_water",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimless/dimTime,0)
);

if (eventIsPresent_water)
{
    event_water.setTimeScheme(h.name(), mesh);
    event_water.setFieldDimensions(sourceTerm_water.dimensions());
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

