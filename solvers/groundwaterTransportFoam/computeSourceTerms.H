//- compute water source term
if (eventIsPresent_water)
{
    sourceTerm_water.primitiveFieldRef() = 0;
    if (isBackwardScheme)
    {
        scalar coefft0_00 = runTime.deltaTValue()/(runTime.deltaTValue() + runTime.deltaT0Value());
        scalar coefftn_0 = 1 + coefft0_00;

        forAll(event_water.coordinates(),pointi)
        {
            sourceTerm_water[idCoordinates_water[pointi]] += (coefftn_0*event_water.currentValue(pointi) - coefft0_00*event_water.oldValue(pointi))/mesh.V()[idCoordinates_water[pointi]];
        }
    }
    else
    {
        forAll(event_water.coordinates(),pointi)
        {
            sourceTerm_water[idCoordinates_water[pointi]] += event_water.currentValue(pointi)/mesh.V()[idCoordinates_water[pointi]];
        }
    }
}

//- compute tracer source term
if (eventIsPresent_tracer)
{
    sourceTerm_tracer.primitiveFieldRef() = 0;
    if (isBackwardScheme)
    {
        scalar coefft0_00 = runTime.deltaTValue()/(runTime.deltaTValue() + runTime.deltaT0Value());
        scalar coefftn_0 = 1 + coefft0_00;

        forAll(event_tracer.coordinates(),pointi)
        {
            sourceTerm_tracer[idCoordinates_tracer[pointi]] += (coefftn_0*event_tracer.currentValue(pointi) - coefft0_00*event_tracer.oldValue(pointi))/mesh.V()[idCoordinates_tracer[pointi]];
        }
    }
    else
    {
        forAll(event_tracer.coordinates(),pointi)
        {
            sourceTerm_tracer[idCoordinates_tracer[pointi]] += event_tracer.currentValue(pointi)/mesh.V()[idCoordinates_tracer[pointi]];
        }
    }
}
