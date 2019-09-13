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

labelList idCoordinates_water(event_water.ncoordinates());
label currentEvent_water = 0;
if (eventIsPresent_water)
{
    sourceTerm_water = dimensionedScalar("",dimless/dimTime,0);
    forAll(event_water.coordinates(),pointi)
    {
        idCoordinates_water[pointi] = mesh.findNearestCell(event_water.coordinates()[pointi]);
        sourceTerm_water[idCoordinates_water[pointi]] = event_water.datas()[currentEvent_water][pointi]/mesh.V()[idCoordinates_water[pointi]];
    }
}

//- Tracer source event reading
word eventFileName_tracer = transportProperties.lookupOrDefault<word>("eventFileTracerSource","");
bool eventIsPresent_tracer = transportProperties.found("eventFileTracerSource");
sourceEventFile event_tracer(eventFileName_tracer);

volScalarField sourceTerm_tracer
(
    IOobject
    (
        "sourceTerm_tracer",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimMass/(dimVol*dimTime),0)
);

labelList idCoordinates_tracer(event_tracer.ncoordinates());
label currentEvent_tracer = 0;
if (eventIsPresent_tracer)
{
    sourceTerm_tracer = dimensionedScalar("",dimMass/(dimVol*dimTime),0);
    forAll(event_tracer.coordinates(),pointi)
    {
        idCoordinates_tracer[pointi] = mesh.findNearestCell(event_tracer.coordinates()[pointi]);
        sourceTerm_tracer[idCoordinates_tracer[pointi]] = event_tracer.datas()[currentEvent_tracer][pointi]/mesh.V()[idCoordinates_tracer[pointi]];
    }
}

//- Output event reading
bool outputEventIsPresent = runTime.controlDict().found("eventFileOutput");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("eventFileOutput","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent) outputEvent.update(runTime.startTime().value());

