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
    Read the control parameters used by setDeltaT

\*---------------------------------------------------------------------------*/

bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

scalar truncationError_global =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError",0.01);

scalar truncationError_C =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError_C",truncationError_global);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

PtrList<timestepManager> dtManagerC(0);
forAll(composition.Y(), speciei)
{
    dtManagerC.append(new timestepManager(runTime, composition.Y()[speciei], truncationError_C));
}

if (adjustTimeStep)
{
    Info << nl << "Timestepping is based on time-scheme truncation error with :"
        << nl << "{"
        << nl << "    truncationError = " << truncationError_C
        << nl << "}"
        << endl;
 }

const bool eventTimeTracking =
    runTime.controlDict().lookupOrDefault("eventTimeTracking", false);
// ************************************************************************* //
