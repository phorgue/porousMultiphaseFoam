{
    //- Set h equal to old-time value
    h = h.oldTime();
    porousModel->rewindTime();
    //- Rewind time
    runTime.setTime(runTime.timeOutputValue()-runTime.deltaTValue(),runTime.timeIndex());
    //- Recompute time step
    iterPicard += 1;
    iterNewton += 1;
    #include "setDeltaT.H"
    //- Update new time
    runTime.setTime(runTime.timeOutputValue()+runTime.deltaTValue(),runTime.timeIndex());
    //- Update properties
    #include "updateProperties.H"
}
