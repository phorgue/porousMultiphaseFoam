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
    //- in case of non convergence of the Picard or Newton  algorithm
    if (iterPicard == maxIterPicard+1 || iterNewton == maxIterNewton+1)
    {
        scalar deltaTFact = dTFactDecrease;
        runTime.setDeltaT
            (
                min(
                    deltaTFact*runTime.deltaTValue(),
                    maxDeltaT
                )
            );
    }
    else
    {
            scalar dtForh = Foam::pow(2*truncationError_h*(hmax+VSMALL)/(dh2dT2max+VSMALL),1./2.);
            scalar newDeltaT = min(dtForh, 1.2*runTime.deltaTValue());
            runTime.setDeltaT
                (
                    min(
                        newDeltaT,
                        maxDeltaT
                    )
                );
    }

    scalar timeOfNextEvent = GREAT;
    if (eventTimeTracking)
    {
        if (sourceEventIsPresent) timeOfNextEvent = min(timeOfNextEvent,sourceEvent.currentEventEndTime());
        forAll(patchEventList,patchEventi) timeOfNextEvent = min(timeOfNextEvent,patchEventList[patchEventi]->currentEventEndTime());
    }

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
        if (eventTimeTracking)
        {
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
        }
        runTime.setDeltaTNoAdjust(min(runTime.deltaTValue(),timeToCloseEvent));
    }

    Info<< "deltaT = " <<  runTime.deltaTValue() << " ; dtheta max = " << dtheta << endl;

}

// ************************************************************************* //
