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

Global
    setDeltaT.H

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

    This routine also adjust timestep to ensure that event times (for patch
    or source term) is exactly reached by the solver.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    scalar dtForC = GREAT;
    forAll(composition.Y(), speciei)
    {
        dtForC = min
            (
                dtForC,
                dtManager[speciei].computeTimestep()
            );
    }

    scalar newDeltaT = min(dtForC, 1.2*runTime.deltaTValue());
    runTime.setDeltaT
        (
            min
            (
                newDeltaT,
                maxDeltaT
            )
        );

    //-Adjust time step to explicitly compute (source/tracer) event time
    if (eventTimeTracking)
    {
        scalar timeOfNextEvent = GREAT;
        forAll(sourceEventList,sourceEventi) timeOfNextEvent = min(timeOfNextEvent,sourceEventList[sourceEventi]->currentEventEndTime());
        forAll(patchEventList,patchEventi) timeOfNextEvent = min(timeOfNextEvent,patchEventList[patchEventi]->currentEventEndTime());

        scalar timeToNextEvent = timeOfNextEvent-runTime.timeOutputValue();
        scalar nSteps =  timeToNextEvent/runTime.deltaTValue();
        if ((nSteps < labelMax) && (nSteps != 0.0))
        {
            const label nStepsToNextEvent = label(max(nSteps, 1) + 0.99);
            runTime.setDeltaT(timeToNextEvent/nStepsToNextEvent,false);
        }

        //- To handle close event times (inferior to current timestep)
        if (nSteps == 0)
        {
            scalar timeToCloseEvent = GREAT;
            forAll(sourceEventList,sourceEventi)
            {
                if (sourceEventList[sourceEventi]->currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,sourceEventList[sourceEventi]->currentEventEndTime()-runTime.timeOutputValue());
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
