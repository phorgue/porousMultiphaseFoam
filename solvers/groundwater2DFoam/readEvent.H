word eventFileName = transportProperties.lookupOrDefault<word>("eventFile","");
bool eventIsPresent = transportProperties.found("eventFile");
eventFile event(eventFileName);
scalar zScale(mesh.bounds().max().z()-mesh.bounds().min().z());

//- determining cell id corresponding to datas
labelList idCoordinates(event.ncoordinates());
label currentEvent = 0;
if (eventIsPresent)
{
    infiltration = dimensionedScalar("",dimLength/dimTime,0);
    forAll(event.coordinates(),pointi)
    {
        idCoordinates[pointi] = mesh.findNearestCell(event.coordinates()[pointi]);
        infiltration[idCoordinates[pointi]] = event.datas()[currentEvent][pointi]*zScale/mesh.V()[idCoordinates[pointi]];
    }
}
