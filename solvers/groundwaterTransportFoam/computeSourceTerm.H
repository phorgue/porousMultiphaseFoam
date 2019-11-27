//- compute water source term
if (eventIsPresent_water)
{
    event_water.updateValue(runTime);
    sourceTerm_water = event_water.dtValuesAsField();
}

//tracer source terms are updated by the mixture class
