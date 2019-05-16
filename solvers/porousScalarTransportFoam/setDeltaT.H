/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    //- C-variation control
    scalar dtForC =  max(dCmax*gMax(C),dCresidual)/(dCdTmax+SMALL);
    scalar newDeltaT = min(dtForC, 1.2*runTime.deltaTValue());
    runTime.setDeltaT
        (
            min
            (
                newDeltaT,
                maxDeltaT
            )
        );

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

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;

}

// ************************************************************************* //
