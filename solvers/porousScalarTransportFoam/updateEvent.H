if (eventIsPresent)
{
    if (currentEvent < event.ndates()-1)
    {
        if (runTime.timeOutputValue() > event.dates()[currentEvent+1])
        {
            currentEvent++;
            sourceTerm = dimensionedScalar("",dimMass/(dimVol*dimTime),0);
            forAll(event.coordinates(),pointi)
            {
                sourceTerm[idCoordinates[pointi]] = event.datas()[currentEvent][pointi]/mesh.V()[idCoordinates[pointi]];
            }
        }
    }
}
