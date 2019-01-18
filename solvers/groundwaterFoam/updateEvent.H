if (eventIsPresent)
{
    if (currentEvent < event.ndates()-1)
    {
        if (runTime.timeOutputValue() > event.dates()[currentEvent+1])
        {
            currentEvent++;
            sourceTerm = dimensionedScalar("",dimLength/dimTime,0);
            forAll(event.coordinates(),pointi)
            {
                Info << idCoordinates[pointi] << endl;
                sourceTerm[idCoordinates[pointi]] = event.datas()[currentEvent][pointi];
            }
        }
    }
}
