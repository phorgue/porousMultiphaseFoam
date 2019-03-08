//- write using event file or controlDict
if (outputEventIsPresent)
{
    if (outputEvent.currentEventStartTime() == runTime.timeOutputValue())
    {
        runTime.writeNow();
    }
}
else
{
    runTime.write();
}
