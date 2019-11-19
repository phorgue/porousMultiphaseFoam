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

//tracer source terms are updated automatically by the mixture class
