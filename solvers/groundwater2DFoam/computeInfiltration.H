if (infiltrationEventIsPresent)
{
    infiltrationEvent.updateValue(runTime);
    if (infiltrationEvent.isUniform()) infiltration.primitiveFieldRef() = infiltrationEvent.dtValue(0);
    else infiltration.primitiveFieldRef() = infiltrationEvent.dtValues();
}

if (waterSourceEventIsPresent)
{
    waterSourceEvent.updateValue(runTime);
    waterSourceTerm = waterSourceEvent.dtValuesAsField();
    Info << nl << "waterSourceTerm min " << min(waterSourceTerm).value() << ", max = " << max(waterSourceTerm).value() << endl;
}
