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
    Timestep for groundwaterFoam solver

\*---------------------------------------------------------------------------*/

if(adjustTimeStep)
{
    //- in case of non convergence of the Picard or Newton  algorithm
    if (iterPicard == maxIterPicard+1 || iterNewton == maxIterNewton+1)
    {
        runTime.setDeltaT(dTFactDecrease*runTime.deltaTValue());
    }
    else if (truncation_mode)
    {
        scalar dtForh = dtManager.computeTimestep();
        scalar newDeltaT = min(dtForh, 1.2*runTime.deltaTValue());
        runTime.setDeltaT(min(newDeltaT, maxDeltaT));
    }
    else
    {
        if (iterPicard < nIterPicard)
        {
            nIterIncrease++;
            if (nIterIncrease == 5)
            {
                runTime.setDeltaT(dTFactIncrease*runTime.deltaTValue());
                nIterIncrease = 0;
            }

        }
        else
        {
            nIterIncrease = 0;
        }
    }

    //-Adjust time step to explicitly compute (source/patch) event time
    if (eventTimeTracking)
    {
        scalar timeOfNextEvent = GREAT;
        if (sourceEventIsPresent) timeOfNextEvent = min(timeOfNextEvent,sourceEvent.currentEventEndTime());
        forAll(patchEventList,patchEventi) timeOfNextEvent = min(timeOfNextEvent,patchEventList[patchEventi]->currentEventEndTime());

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
            if (sourceEventIsPresent)
            {
                if (sourceEvent.currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,sourceEvent.currentEventEndTime()-runTime.timeOutputValue());
                }
            }
            forAll(patchEventList,patchEventi)
            {
                if (patchEventList[patchEventi]->currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,patchEventList[patchEventi]->currentEventEndTime()-runTime.timeOutputValue());
                }
            }
        runTime.setDeltaT(min(runTime.deltaTValue(),timeToCloseEvent),false);
        }
    }

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;

}

// ************************************************************************* //
