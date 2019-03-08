//- write using event file or controlDict
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() == runTime.timeOutputValue())
    {
        runTime.writeNow();
    }
}
else
{
    runTime.write();
}
