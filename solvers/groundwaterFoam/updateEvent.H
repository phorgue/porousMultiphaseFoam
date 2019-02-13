if (eventIsPresent)
{
    if (currentEvent < (event.ndates()-1))
    {
        sourceTerm = dimensionedScalar("",dimless/dimTime,0);
        if (runTime.timeOutputValue() > event.dates()[currentEvent+1])
        {
            scalar eventFactor =  (runTime.timeOutputValue()- event.dates()[currentEvent+1]) / runTime.deltaTValue();
        forAll(event.coordinates(),pointi)
        {
            sourceTerm[idCoordinates[pointi]] = (eventFactor * event.datas()[currentEvent+1][pointi]
            + (1-eventFactor) * event.datas()[currentEvent][pointi]) / mesh.V()[idCoordinates[pointi]];
        }
        currentEvent++;
    }
    else
    {
        forAll(event.coordinates(),pointi)
        {
            sourceTerm[idCoordinates[pointi]] = event.datas()[currentEvent][pointi]/mesh.V()[idCoordinates[pointi]];
        }
    }
}
}
