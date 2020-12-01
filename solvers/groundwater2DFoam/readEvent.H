
////////////////////////////////////////////////////////
////////////////// INFILTRATION EVENT //////////////////
////////////////////////////////////////////////////////

bool infiltrationEventIsPresent = transportProperties.found("eventFileInfiltration");
word infiltrationEventFileName = transportProperties.lookupOrDefault<word>("eventFileInfiltration","");
infiltrationEventFile infiltrationEvent(infiltrationEventFileName);
if (infiltrationEventIsPresent)
{
    infiltrationEvent.updateIndex(runTime.startTime().value());
    infiltrationEvent.updateValue(runTime);
    infiltrationEvent.storeOldValues();
}

scalar eventTimeStep(runTime.controlDict().lookupOrDefault<scalar>("eventTimeStep",0));
if (eventTimeStep > 0)
{
    if (infiltrationEventIsPresent) infiltrationEvent.addIntermediateTimeSteps(eventTimeStep);
}

//- Updating infiltration value using event file
if (infiltrationEventIsPresent)
{
    infiltrationEvent.setTimeScheme(potential.name(), mesh);

    if (infiltrationEvent.isUniform())
    {
        infiltration.primitiveFieldRef() = infiltrationEvent.currentValue(0);
    }
    else
    {
        if (infiltrationEvent.currentValues().size() != mesh.C().size())
        {
            Warning() << "Number of infiltration values does not correspond to the mesh size, first value is used as uniform infiltration" << endl;
            infiltrationEvent.setInfiltrationAsUniform();
            infiltration.primitiveFieldRef() = infiltrationEvent.currentValue(0);
        }
        else
        {
            infiltration.primitiveFieldRef() = infiltrationEvent.currentValues();
        }

    }
}

//////////////////////////////////////////////////////
//////////////////// SOURCE EVENT ////////////////////
//////////////////////////////////////////////////////

bool waterSourceEventIsPresent = transportProperties.found("eventFileWaterSource");
word waterSourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileWaterSource","");
sourceEventFile waterSourceEvent(waterSourceEventFileName);
if (waterSourceEventIsPresent)
{
    waterSourceEvent.setTimeScheme(potential.name(), mesh);
    waterSourceEvent.setFieldDimensions(waterSourceTerm.dimensions());
    waterSourceEvent.updateIndex(runTime.startTime().value());
    waterSourceEvent.updateValue(runTime);
    waterSourceEvent.storeOldValues();
}

//////////////////////////////////////////////////////
//////////////////// OUTPUT EVENT ////////////////////
//////////////////////////////////////////////////////

bool outputEventIsPresent = runTime.controlDict().found("eventFileOutput");
word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("eventFileOutput","");
outputEventFile outputEvent(outputEventFileName);
if (outputEventIsPresent) outputEvent.updateIndex(runTime.startTime().value());

bool CSVoutput = runTime.controlDict().lookupOrDefault<bool>("CSVoutput",true);
OFstream waterMassBalanceCSV("waterMassBalance.csv");
if (CSVoutput)
{
    waterMassBalanceCSV << "#Time flux(Infiltration)";
    if (waterSourceEventIsPresent) waterMassBalanceCSV << " flux(waterSource)";
    if (fixedPotentialIDList.size() > 0) waterMassBalanceCSV << " flux(fixedPoints)";
    if (seepageActive) waterMassBalanceCSV << " flux(seepage)";
    forAll(mesh.boundaryMesh(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            waterMassBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
        }
    }
    waterMassBalanceCSV << endl;
}
