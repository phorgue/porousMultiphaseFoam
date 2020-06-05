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

//- Global parameters
scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
const word timeStepControl =
    runTime.controlDict().lookupOrDefault<word>("timeStepControl", "Picard");
const bool eventTimeTracking =
     runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

//- if dtFromPicard used
const label nIterStability = runTime.controlDict().lookupOrDefault<label>("nIterStability",5);
const label minIterPicard = runTime.controlDict().lookupOrDefault<label>("minIter",3);
label iterStability = 0;
const scalar dTFactIncrease = runTime.controlDict().lookupOrDefault<scalar>("dTFactIncrease",1.2);
const scalar dTFactDecrease = runTime.controlDict().lookupOrDefault<scalar>("dTFactDecrease",0.8);

//- if dtFromtheta used
const scalar dthetamax = runTime.controlDict().lookupOrDefault("dthetamax",0.01);
scalar dthetadTmax = dthetamax/runTime.deltaTValue();
scalar dtheta = 0;
scalar dtheta_avg = 0;

if (adjustTimeStep)
{
    if (timeStepControl == "Picard")
    {
        Info << "Timestepping is based on Picard's iteration with :"
            << nl << "{"
            << nl << "    dTFactIncrease = " << dTFactIncrease
            << nl << "    dTFactDecrease = " << dTFactDecrease
            << nl << "    minIterPicard  = " <<  minIterPicard
            << nl << "    nIterStability = " << nIterStability
            << nl << "}"
            << endl;
    }
    else if (timeStepControl == "dthetamax")
    {
        Info << "Timestepping is based on theta variation with :"
            << nl << "{"
            << nl << "    dthetamax = " << dthetamax
            << nl << "    dTFactDecrease = " << dTFactDecrease
            << nl << "}"
            << endl;
    }
    else
    {
        FatalErrorIn("readPicardControls.H") << " timeStepControl is not defined => " << timeStepControl << abort(FatalError);
    }       
}
else
{
    Info << "Timestep is fixed";
}

// ************************************************************************* //
