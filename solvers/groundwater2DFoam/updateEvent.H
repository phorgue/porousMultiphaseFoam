if (eventIsPresent)
{
    if (currentEvent < event.ndates()-1)
    {
        if (runTime.timeOutputValue() > event.dates()[currentEvent+1])
        {
            currentEvent++;
            infiltration = dimensionedScalar("",dimLength/dimTime,0);
            forAll(event.coordinates(),pointi)
            {
                infiltration[idCoordinates[pointi]] = event.datas()[currentEvent][pointi]*zScale/mesh.V()[idCoordinates[pointi]];
            }
        }
    }
}
