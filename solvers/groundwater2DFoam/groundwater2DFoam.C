/*---------------------------------------------------------------------------* \
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    groundwater2DFoam

Description
    Transient solver for free-surface flow in porous media

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "harmonic.H"
#include "fixedValueFvPatchField.H"
#include "DEMfile.H"
#include "infiltrationEventFile.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "timestepManager.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addBoolOption("steady", "to run steady flow simulation");

    Foam::argList args(argc, argv);
    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }
    #include "../headerPMF.H"
    bool steady = args.found("steady");

    Info << "Create time\n" << Foam::endl;
    Time runTime(Time::controlDictName, args);

    #include "createMesh.H"
    #include "createFields.H"
    #include "readFixedPoints.H"
    #include "readTimeControls.H"

    autoPtr<infiltrationEventFile> infiltrationEvent = infiltrationEventFile::New("infiltrationEventFile", transportProperties);
    infiltrationEvent->init(runTime, potential.name(), mesh, infiltration);
    autoPtr<sourceEventFile> sourceEvent = sourceEventFile::New("sourceEventFileWater", transportProperties);
    sourceEvent->init(runTime, potential.name(), mesh, waterSourceTerm.dimensions());
    autoPtr<outputEventFile> outputEvent = outputEventFile::New(runTime, zScale);
    outputEvent->addField(hwater, phi, eps, "m3", true);
    outputEvent->addField(potential, phi, "m", false);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        if (!steady)
        {
            if (infiltrationEvent->isPresent()) infiltrationEvent->updateIndex(runTime.timeOutputValue());
            if (sourceEvent->isPresent()) sourceEvent->updateIndex(runTime.timeOutputValue());
            #include "setDeltaT.H"
        }

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        //- Update infiltration term
        if (!steady)
        {
            if (infiltrationEvent->isPresent()) infiltrationEvent->updateInfiltration(runTime, infiltration.primitiveFieldRef());

           if (sourceEvent->isPresent()) {
                sourceEvent->updateValue(runTime);
                waterSourceTerm = sourceEvent->dtValuesAsField();
            }
        }

        //- Solve potential equation
        #include "potentialEqn.H"

        //- Residual computation
        if (steady)
        {
            if (maxResidual < residualPotential)
            {
                runTime.writeAndEnd();
            }
            else
            {
                runTime.write();
            }
        }
        else
        {
            outputEvent->write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    if (cumulativeWaterAdded > 0) Info << "Cumulated water added = " << cumulativeWaterAdded << " m3, equivalent height = " << cumulativeWaterAdded*zScale/gSum(mesh.V()) << " m" << nl << endl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
