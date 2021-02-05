//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() <= runTime.timeOutputValue())
    {
        //- compute and write time interpolated fields
        volScalarField p_interpolated = outputEvent.timeInterpolate(p, runTime);
        volScalarField Sb_interpolated = outputEvent.timeInterpolate(Sb, runTime);
        volVectorField Ua_interpolated = outputEvent.timeInterpolate(Ua, runTime);
        volVectorField Ub_interpolated = outputEvent.timeInterpolate(Ub, runTime);
        surfaceScalarField phia_interpolated = outputEvent.timeInterpolate(phia, runTime);
        surfaceScalarField phib_interpolated = outputEvent.timeInterpolate(phib, runTime);

        outputEvent.updateIndex(runTime.timeOutputValue());
    }
}
else
{
    //- write fields using usual openfoam rules
    runTime.write();
}
