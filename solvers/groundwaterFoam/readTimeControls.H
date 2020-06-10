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

Global
    readTimeControls

Description
    Read the control parameters used by setDeltaT for groundwaterFoam solver
\*---------------------------------------------------------------------------*/

const scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

const scalar truncationError_global =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError",0.01);

const scalar truncationError_h =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError_h",truncationError_global);

const scalar dTFactDecrease =
    runTime.controlDict().lookupOrDefault<scalar>("dTFactDecrease",0.8);

Foam::fv::EulerD2dt2Scheme<scalar> d2dt2Operator(mesh);
scalar dtheta = 0;
scalar hmax = 0;
scalar dh2dT2max = 0;

if (adjustTimeStep)
{
    //- Derivative initialisation to keep 1st user-defined time step
    if (adjustTimeStep)
    {
        hmax = gMax(h);
        dh2dT2max = 2*truncationError_h*(hmax+VSMALL)/Foam::pow(runTime.deltaTValue(),2);
        Info << "Timestepping is based on time-scheme truncation error with :"
            << nl << "{"
            << nl << "    truncationError = " << truncationError_h
            << nl << "    dTFactDecrease  = " << dTFactDecrease
            << nl << "}"
            << endl;
    }
    else
    {
        Info << "Timestep is fixed";
    }
}

const bool eventTimeTracking =
     runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

// ************************************************************************* //
