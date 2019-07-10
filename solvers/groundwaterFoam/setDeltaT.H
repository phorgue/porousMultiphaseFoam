/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

    //- Timestep using number of Picard's iterations
    if (timeStepControl == "Picard")
    {
        scalar deltaTFact=1;
        if(iterPicard == maxIterPicard)
        {
            deltaTFact = dTFactDecrease;
            iterStability=0;
        }
        else if (iterPicard <= minIterPicard)
        {
            iterStability++;
        }

        if (iterStability == nIterStability)
        {
            deltaTFact = dTFactIncrease;
            iterStability = 0;
        }

        runTime.setDeltaT
            (
                min(
                    deltaTFact*runTime.deltaTValue(),
                    maxDeltaT
                )
            );
    }

    //- Timestep using h variation
    else if(timeStepControl == "dthetamax")
    {
        scalar dtForTheta = dthetamax / (dthetadTmax+SMALL);

        runTime.setDeltaT
            (
                min(
                    dtForTheta,
                    min(
                        1.2*runTime.deltaTValue(),
                        maxDeltaT
                    )
                )
            );
    }

    scalar timeOfNextEvent = GREAT;
    if (sourceEventIsPresent) timeOfNextEvent = min(timeOfNextEvent,sourceEvent.currentEventEndTime());
    if (patchEventIsPresent) timeOfNextEvent = min(timeOfNextEvent,patchEvent.currentEventEndTime());
    if (outputEventIsPresent) timeOfNextEvent = min(timeOfNextEvent,outputEvent.currentEventEndTime());

    scalar timeToNextEvent = timeOfNextEvent-runTime.timeOutputValue();
    scalar nSteps =  timeToNextEvent/runTime.deltaTValue();

    if ((nSteps < labelMax) && (nSteps != 0))
    {
        const label nStepsToNextEvent = label(max(nSteps, 1) + 0.99);
        runTime.setDeltaTNoAdjust(timeToNextEvent/nStepsToNextEvent);
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
        if (patchEventIsPresent)
        {
            if (patchEvent.currentEventEndTime() != runTime.timeOutputValue())
            {
                timeToCloseEvent = min(timeToCloseEvent,patchEvent.currentEventEndTime()-runTime.timeOutputValue());
            }
        }
        if (outputEventIsPresent)
        {
            if (outputEvent.currentEventEndTime() != runTime.timeOutputValue())
            {
                timeToCloseEvent = min(timeToCloseEvent,outputEvent.currentEventEndTime()-runTime.timeOutputValue());
            }
        }
        runTime.setDeltaTNoAdjust(min(runTime.deltaTValue(),timeToCloseEvent));
    }

    Info<< "deltaT = " <<  runTime.deltaTValue() << " ; dtheta max = " << dtheta << " ; dtheta avg = " << dtheta_avg << endl;

}

// ************************************************************************* //
