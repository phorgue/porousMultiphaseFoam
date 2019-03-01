bool sourceEventIsPresent = transportProperties.found("eventFileSourceMassFlowRate");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileSourceMassFlowRate","");
sourceEventFile sourceEvent(sourceEventFileName);
if (sourceEventIsPresent) sourceEvent.update(runTime.startTime().value());

bool patchEventIsPresent =  transportProperties.found("eventFilePatchMassFlowRate");
word patchEventFileName = transportProperties.lookupOrDefault<word>("eventFilePatchMassFlowRate","");
patchEventFile patchEvent(patchEventFileName,false);
if (patchEventIsPresent) patchEvent.update(runTime.startTime().value());

bool outputEventIsPresent = runTime.controlDict().found("eventFileOutput");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("eventFileOutput","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent) outputEvent.update(runTime.startTime().value());


scalar eventTimeStep(runTime.controlDict().lookupOrDefault<scalar>("eventTimeStep",0));
if (eventTimeStep > 0)
{
    if (patchEventIsPresent) patchEvent.addIntermediateTimeSteps(eventTimeStep);
    if (sourceEventIsPresent) sourceEvent.addIntermediateTimeSteps(eventTimeStep);
}

//- Computing source term
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

//- determining cell id corresponding to datas
labelList idCoordinates(sourceEvent.ncoordinates());
if (sourceEventIsPresent)
{
    sourceTerm = dimensionedScalar("",dimMass/(dimVol*dimTime),0);
    forAll(sourceEvent.coordinates(),pointi)
    {
        idCoordinates[pointi] = mesh.findNearestCell(sourceEvent.coordinates()[pointi]);
        sourceTerm[idCoordinates[pointi]] = sourceEvent.currentValue(pointi)/mesh.V()[idCoordinates[pointi]];
    }
}
