bool sourceEventIsPresent = transportProperties.found("eventFileSourceMassFlowRate");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileSourceMassFlowRate","");

bool patchEventIsPresent =  transportProperties.found("eventFilePatchMassFlowRate");
word patchEventFileName = transportProperties.lookupOrDefault<word>("eventFilePatchMassFlowRate","");

sourceEventFile sourceEvent(sourceEventFileName);
patchEventFile patchEvent(patchEventFileName,false);

volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimMass/(dimVol*dimTime),0)
);

//- Initialising event iterators
label currentSourceEvent = 0;
label currentPatchEvent = 0;
label currentEventTimeStep = 0;
label currentEventWrite = 1;

//- determining cell id corresponding to datas
labelList idCoordinates(sourceEvent.ncoordinates());
if (sourceEventIsPresent)
{
    sourceTerm = dimensionedScalar("",dimMass/(dimVol*dimTime),0);
    forAll(sourceEvent.coordinates(),pointi)
    {
        idCoordinates[pointi] = mesh.findNearestCell(sourceEvent.coordinates()[pointi]);
        sourceTerm[idCoordinates[pointi]] = sourceEvent.datas()[currentSourceEvent][pointi]/mesh.V()[idCoordinates[pointi]];
    }
}

//- Constructing list of all event times
scalarList allEventDates(sourceEvent.ndates()+patchEvent.ndates());

if (patchEventIsPresent && sourceEventIsPresent)
{
    label iterPatch = 0;
    label iterSource = 0;
    label iterGlobal = 0;
    allEventDates = 0;
    allEventDates[0] = runTime.startTime().value();
    while ((iterPatch < patchEvent.ndates()) && (iterSource < sourceEvent.ndates()))
    {
        if (patchEvent.dates()[iterPatch] <= sourceEvent.dates()[iterSource])
        {
            if(patchEvent.dates()[iterPatch] > allEventDates[iterGlobal])
            {
                iterGlobal++;
                allEventDates[iterGlobal] = patchEvent.dates()[iterPatch];
            }
            iterPatch++;
        }
        else
        {
            if(sourceEvent.dates()[iterSource] > allEventDates[iterGlobal])
            {
                iterGlobal++;
                allEventDates[iterGlobal] = sourceEvent.dates()[iterSource];
            }
            iterSource++;

        }
    }

    allEventDates.resize(iterGlobal);

}
else if (patchEventIsPresent)
{
    allEventDates = patchEvent.dates();
}
else if (sourceEventIsPresent)
{
    allEventDates = sourceEvent.dates();
}

Info << nl << "There are " << allEventDates.size() << " events for this simulation" << endl;

//- Adding intermediate event to force small time step.
scalar eventTimeStep(runTime.controlDict().lookupOrDefault<scalar>("eventTimeStep",1));
scalarList allEventDatesForTimeStep = allEventDates;

if (allEventDates.size() > 2)
{
    label nbEvents = allEventDates.size();
    allEventDatesForTimeStep.resize((nbEvents-2)*3+2);
    allEventDatesForTimeStep[0] = allEventDates[0];
    for(label eventi=1;eventi<nbEvents-1;eventi++)
    {
        allEventDatesForTimeStep[eventi*3-2] = allEventDates[eventi]-eventTimeStep;
        allEventDatesForTimeStep[eventi*3-1] = allEventDates[eventi];
        allEventDatesForTimeStep[eventi*3] = allEventDates[eventi]+eventTimeStep;
    }
    allEventDatesForTimeStep[(nbEvents-2)*3+1] =  allEventDates[nbEvents-1];
}
