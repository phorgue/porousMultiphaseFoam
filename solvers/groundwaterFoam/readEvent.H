word eventFileName = transportProperties.lookupOrDefault<word>("eventFile","");
bool eventIsPresent = transportProperties.found("eventFile");
eventFile event(eventFileName);

volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimVol/dimTime,0)
);

//- determining cell id corresponding to datas
labelList idCoordinates(event.ncoordinates());
label currentEvent = 0;
if (eventIsPresent)
{
    sourceTerm = dimensionedScalar("",dimVol/dimTime,0);
    forAll(event.coordinates(),pointi)
    {
        idCoordinates[pointi] = mesh.findNearestCell(event.coordinates()[pointi]);
        sourceTerm[idCoordinates[pointi]] = event.datas()[currentEvent][pointi];
    }
}
