/*---------------------------------------------------------------------------*\
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

Description
    Timestep for groundwater2DFoam solver

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    scalar newDeltaT = min(dtManager.computeTimestep(),1.2*runTime.deltaTValue());

    runTime.setDeltaT
        (
            min
            (
                newDeltaT,
                maxDeltaT
            )
        );

    //-Adjust time step to explicitly compute infiltration event time
    if (eventTimeTracking)
    {
        scalar timeOfNextEvent = GREAT;
        if (infiltrationEventIsPresent) timeOfNextEvent = min(timeOfNextEvent,infiltrationEvent.currentEventEndTime());

        scalar timeToNextEvent = timeOfNextEvent-runTime.timeOutputValue();
        scalar nSteps =  timeToNextEvent/runTime.deltaTValue();
        if ((nSteps < labelMax) && (nSteps != 0))
        {
            const label nStepsToNextEvent = label(max(nSteps, 1) + 0.99);
            runTime.setDeltaT(timeToNextEvent/nStepsToNextEvent,false);
        }

        //- To handle close event times (inferior to current timestep)
        if (nSteps == 0)
        {
            scalar timeToCloseEvent = GREAT;
            if (infiltrationEventIsPresent)
            {
                if (infiltrationEvent.currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,infiltrationEvent.currentEventEndTime()-runTime.timeOutputValue());
                }
            }
            runTime.setDeltaT(min(runTime.deltaTValue(),timeToCloseEvent),false);
        }
    }

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
           
}    


// ************************************************************************* //
