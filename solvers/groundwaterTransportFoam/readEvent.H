word eventFileName = transportProperties.lookupOrDefault<word>("eventFileFlowRate","");
bool eventIsPresent = transportProperties.found("eventFileFlowRate");
sourceEventFile event(eventFileName);

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
    dimensionedScalar("",dimless/dimTime,0)
);

//- determining cell id corresponding to datas
labelList idCoordinates(event.ncoordinates());
label currentEvent = 0;
if (eventIsPresent)
{
    sourceTerm = dimensionedScalar("",dimless/dimTime,0);
    forAll(event.coordinates(),pointi)
    {
        idCoordinates[pointi] = mesh.findNearestCell(event.coordinates()[pointi]);
        sourceTerm[idCoordinates[pointi]] = event.datas()[currentEvent][pointi]/mesh.V()[idCoordinates[pointi]];
    }
}
