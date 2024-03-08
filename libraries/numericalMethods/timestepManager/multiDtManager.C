
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

\*---------------------------------------------------------------------------*/

#include "multiDtManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

multiDtManager::multiDtManager(
    Time& runTime,
    const List<sourceEventFile*>& sourceEventList,
    const List<patchEventFile*>& patchEventList
)
    :
    runTime_(runTime),
    dtManagerT_(),
    sourceEventList_(sourceEventList),
    patchEventList_(patchEventList),
    infiltrationEventList_(0)
{
    adjustTimeStep_ =
        runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);
    maxDeltaT_ =
        runTime_.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
    eventTimeTracking_ =
        runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

     Info << nl << "General time-stepping "
     << nl << "{"
     << nl << "    adjustTimeStep is " << adjustTimeStep_
     << nl << "    maxDeltaT =  " << maxDeltaT_
     << nl << "    eventTimeTracking is " << eventTimeTracking_
     << nl << "}" << endl;

}

multiDtManager::multiDtManager(
    Time& runTime,
    const List<sourceEventFile*>& sourceEventList,
    const List<infiltrationEventFile*>& infiltrationEventList
)
    :
    runTime_(runTime),
    dtManagerT_(),
    sourceEventList_(sourceEventList),
    patchEventList_(0),
    infiltrationEventList_(infiltrationEventList)
{
    adjustTimeStep_ =
        runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);
    maxDeltaT_ =
        runTime_.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
    eventTimeTracking_ =
        runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

     Info << nl << "General time-stepping "
     << nl << "{"
     << nl << "    adjustTimeStep is " << adjustTimeStep_
     << nl << "    maxDeltaT =  " << maxDeltaT_
     << nl << "    eventTimeTracking is " << eventTimeTracking_
     << nl << "}" << endl;

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

multiDtManager::~multiDtManager()
= default;

// * * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * //

void multiDtManager::addField
(
    const volScalarField& field,
    const labelList* fixedCells
)
{
    dtManagerT_.append(new timestepManagerTruncation(runTime_, field, fixedCells));
}

void multiDtManager::addField
(
    const volScalarField& field,
    const labelList& fixedCells
)
{
    const labelList* fixedCellsPtr = &fixedCells;
    dtManagerT_.append(new timestepManagerTruncation(runTime_, field, fixedCellsPtr));
}

    void multiDtManager::addIterativeAlgorithm
(
    const volScalarField& field,
    const word& algoName
)
{
    dtManagerI_.append(new timestepManagerIterative(runTime_, field.mesh().solutionDict(), algoName));
}

void multiDtManager::updateDt()
{
    if (adjustTimeStep_) {
        scalar dt = GREAT;
        forAll(dtManagerT_, fieldi) {
            dt = min
                (
                    dt,
                    dtManagerT_[fieldi].computeTimestep()
                );
        }
        forAll(dtManagerI_, algoi) {
            dt = min
                (
                    dt,
                    dtManagerI_[algoi].computeTimestep()
                );
        }

        dt = min(dt, 1.25 * runTime_.deltaTValue());
        runTime_.setDeltaT(min(dt, maxDeltaT_));
        if (eventTimeTracking_) adjustDeltaTUsingEvent();
        Info << "deltaT = " << runTime_.deltaTValue() << endl;
    }
}

void multiDtManager::updateAllDerivatives()
{
    forAll(dtManagerT_, fieldi) dtManagerT_[fieldi].updateDerivatives();
}


void multiDtManager::adjustDeltaTUsingEvent()
{
    //-Adjust time step to explicitly compute (source/tracer) event time
    scalar timeOfNextEvent = GREAT;
    forAll(sourceEventList_,sourceEventi) timeOfNextEvent = min(timeOfNextEvent,sourceEventList_[sourceEventi]->currentEventEndTime());
    forAll(patchEventList_,patchEventi) timeOfNextEvent = min(timeOfNextEvent,patchEventList_[patchEventi]->currentEventEndTime());
    forAll(infiltrationEventList_,eventi) timeOfNextEvent = min(timeOfNextEvent,patchEventList_[eventi]->currentEventEndTime());

    scalar timeToNextEvent = timeOfNextEvent-runTime_.timeOutputValue();
    scalar nSteps =  timeToNextEvent/runTime_.deltaTValue();
    if ((nSteps < labelMax) && (nSteps != 0.0))
    {
        const auto nStepsToNextEvent = label(max(nSteps, 1) + 0.99);
        runTime_.setDeltaT(timeToNextEvent/nStepsToNextEvent,false);
    }

    //- To handle close event times (inferior to current timestep)
    if (nSteps == 0)
    {
        scalar timeToCloseEvent = GREAT;
        forAll(sourceEventList_,sourceEventi)
        {
            if (sourceEventList_[sourceEventi]->currentEventEndTime() != runTime_.timeOutputValue())
            {
                timeToCloseEvent = min(timeToCloseEvent,sourceEventList_[sourceEventi]->currentEventEndTime()-runTime_.timeOutputValue());
            }
        }
        forAll(patchEventList_,patchEventi)
        {
            if (patchEventList_[patchEventi]->currentEventEndTime() != runTime_.timeOutputValue())
            {
                timeToCloseEvent = min(timeToCloseEvent,patchEventList_[patchEventi]->currentEventEndTime()-runTime_.timeOutputValue());
            }
        }
        runTime_.setDeltaT(min(runTime_.deltaTValue(),timeToCloseEvent),false);
    }
}

} // End namespace Foam

// ************************************************************************* //
