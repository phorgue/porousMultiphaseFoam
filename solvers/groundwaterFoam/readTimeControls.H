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
    readTimeControls

Description
    Read the control parameters used by setDeltaT for groundwaterFoam solver

\*---------------------------------------------------------------------------*/

//-global parameters
const scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

const scalar dTFactDecrease =
    runTime.controlDict().lookupOrDefault<scalar>("dTFactDecrease",0.8);

const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault<bool>("adjustTimeStep", false);

const word timeStepMode =
    runTime.controlDict().lookupOrDefault<word>("timeStepMode", "truncation");

bool truncation_mode = true;
if (timeStepMode == "nIterPicard") truncation_mode = false;
else if (timeStepMode != "truncation")
{
    FatalErrorIn("readTimeControls.H") << "timeStepMode " << timeStepMode << " unknown: should be truncation or nIterPicard"
        << exit(FatalError);
}

//-truncation parameters
const scalar truncationError_global =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError",0);
const scalar truncationError_h =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError_h",truncationError_global);

//-nIterPicard parameters
const label nIterPicard =
    runTime.controlDict().lookupOrDefault<scalar>("nIterPicard",0);
const scalar dTFactIncrease =
    runTime.controlDict().lookupOrDefault<scalar>("dTFactIncrease",1.25);
label nIterIncrease = -1;

if (adjustTimeStep)
{
    if (truncation_mode)
    {
        if (truncationError_h == 0)
        {
            FatalErrorIn("readTimeControls.H") << "truncationError is not defined for timestepping"
        << exit(FatalError);
        }
    }
    else
    {
        if (nIterPicard == 0)
        {
            FatalErrorIn("readTimeControls.H") << "nIterPicard is not defined for timestepping"
        << exit(FatalError);
        }
    }
}


labelList*fixedPotentialIDListPtr  = &fixedPotentialIDList;
timestepManager dtManager(runTime, h, truncationError_h, fixedPotentialIDListPtr);

//- Derivative initialisation to keep 1st user-defined time step
if (adjustTimeStep)
{
    if (truncation_mode)
    {
        Info << nl << "Timestepping is based on time-scheme truncation error with :"
            << nl << "{"
            << nl << "    truncationError = " << truncationError_h
            << nl << "    dTFactDecrease  = " << dTFactDecrease
            << nl << "}"
            << endl;
    }
    else
    {
        Info << nl << "Timestepping is based on number of Picard iterations with :"
            << nl << "{"
            << nl << "    nIterPicard = " << nIterPicard
            << nl << "    dTFactIncrease  = " << dTFactIncrease
            << nl << "    dTFactDecrease  = " << dTFactDecrease
            << nl << "}"
            << endl;
    }
}
else
{
    Info << nl << "Timestep is fixed" << endl;
}

const bool eventTimeTracking =
     runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

// ************************************************************************* //
