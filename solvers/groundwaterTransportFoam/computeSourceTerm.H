//- compute water source term
if (eventIsPresent_water)
{
    sourceTerm_water = event_water.dtValuesAsField();
}

//tracer source terms are updated by the mixture class
