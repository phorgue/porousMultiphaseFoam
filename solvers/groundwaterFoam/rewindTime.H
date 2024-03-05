{
    //- Set h equal to old-time value
    h = h.oldTime();
    porousModel->rewindTime();

    //- Rewind time
    runTime.setTime(runTime.timeOutputValue()-runTime.deltaTValue(),runTime.timeIndex());

    //- Increase iterator for indicate non-convergence
    Picard++;
    Newton++;

    //- Update time step and time output value
    MDTM.updateDt();
    runTime.setTime(runTime.timeOutputValue()+runTime.deltaTValue(),runTime.timeIndex());
    //- Update properties
    #include "updateProperties.H"
}
