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
    else
    {
        scalar dtForh = runTime.deltaTValue();
        //- h-variation control
        if (truncation_mode)
        {
            dtForh = dtManager.computeTimestep();
        }
        else
        {
            if (iterPicard < nIterPicard)
            {
                nIterIncrease++;
                if (nIterIncrease == 5)
                {
                    dtForh = dTFactIncrease*runTime.deltaTValue();
                    nIterIncrease = 0;
                }

            }
            else
            {
                nIterIncrease = 0;
            }
        }

        //- C-variation control
        scalar dtForC = VGREAT;
        forAll(composition.Y(), speciei)
        {
            dtForC = min
                (
                    dtForC,
                    dtManagerC[speciei].computeTimestep()
                );
        }

        //- set minimal delta t
        runTime.setDeltaT
            (
                min(
                    min(dtForh, dtForC),
                    min( 1.2*runTime.deltaTValue(),maxDeltaT)
                )
            );

    }

    //-Adjust time step to explicitly compute (water/tracer) event time 
    scalar timeOfNextEvent = GREAT;
    if (eventTimeTracking)
    {
        if (eventIsPresent_water) timeOfNextEvent = min(timeOfNextEvent,event_water.currentEventEndTime());
        forAll(tracerSourceEventList,tracerSourceEventi) timeOfNextEvent = min(timeOfNextEvent,tracerSourceEventList[tracerSourceEventi]->currentEventEndTime());

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
            if (eventIsPresent_water)
            {
                if (event_water.currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,event_water.currentEventEndTime()-runTime.timeOutputValue());
                }
            }
            forAll(tracerSourceEventList,tracerSourceEventi)
            {
                if (tracerSourceEventList[tracerSourceEventi]->currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,tracerSourceEventList[tracerSourceEventi]->currentEventEndTime()-runTime.timeOutputValue());
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
