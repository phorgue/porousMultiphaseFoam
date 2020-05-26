//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() < runTime.timeOutputValue())
    {
        scalar interpolateFactor = (runTime.timeOutputValue()-outputEvent.currentEventEndTime())/ runTime.deltaTValue();
        scalar timeOutputBackup = runTime.timeOutputValue();
        runTime.setTime(outputEvent.currentEventEndTime(),runTime.timeIndex());
        volScalarField p_interpolated
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                p
            );
        p_interpolated = interpolateFactor*p+(1.0-interpolateFactor)*p.oldTime();
        p_interpolated.write();
        volScalarField Sb_interpolated
        
            (
                IOobject
                (
                    "Sb",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Sb
            );
        Sb_interpolated = interpolateFactor*Sb+(1.0-interpolateFactor)*Sb.oldTime();
        Sb_interpolated.write();
        volVectorField Ua_interpolated
            (
                IOobject
                (
                    "Ua",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Ua
            );
        Ua_interpolated = interpolateFactor*Ua+(1.0-interpolateFactor)*Ua.oldTime();
        Ua_interpolated.write();
        surfaceScalarField phia_interpolated
            (
                IOobject
                (
                    "phia",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phia
            );
        phia_interpolated = interpolateFactor*phia+(1.0-interpolateFactor)*phia.oldTime();
        phia_interpolated.write();

        volVectorField Ub_interpolated
            (
                IOobject
                (
                    "Ub",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Ub
            );
        Ub_interpolated = interpolateFactor*Ub+(1.0-interpolateFactor)*Ub.oldTime();
        Ub_interpolated.write();
        surfaceScalarField phib_interpolated
            (
                IOobject
                (
                    "phib",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phib
            );
        phib_interpolated = interpolateFactor*phib+(1.0-interpolateFactor)*phib.oldTime();
        phib_interpolated.write();

        runTime.setTime(timeOutputBackup,runTime.timeIndex());
        outputEvent.updateIndex(runTime.timeOutputValue());
    }
}
else
{
    //- write fields using usual openfoam rules
    runTime.write();
}
