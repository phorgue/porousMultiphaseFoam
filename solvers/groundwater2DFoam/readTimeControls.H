/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Read the control parameters used by setDeltaT

\*---------------------------------------------------------------------------*/

bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

scalar truncationError =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError",0.01);

word timeScheme = mesh.ddtScheme("ddt(potential)");

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

EulerD3dt3Scheme d3dt3Operator(mesh,runTime.deltaTValue());
Foam::fv::EulerD2dt2Scheme<scalar> d2dt2Operator(mesh);
scalar dpotentialdTmax = 0;
scalar dpotential2dT2max = 0;
scalar dpotential3dT3max = 0;
scalar hwatermax = 0;
//- Derivative initialisation to keep 1st user-defined time step
if (timeScheme == "backward")
{
        hwatermax = gMax(hwater);
        dpotential3dT3max = 3*truncationError*(hwatermax+VSMALL)/Foam::pow(runTime.deltaTValue(),3);
}
else if  (timeScheme == "CrankNicolson")
{
        hwatermax = gMax(hwater);
        dpotential3dT3max = 12*truncationError*(hwatermax+VSMALL)/Foam::pow(runTime.deltaTValue(),3);
}
else if (timeScheme == "Euler")
{
        hwatermax = gMax(hwater);
        dpotential2dT2max = 2*truncationError*(hwatermax+VSMALL)/Foam::pow(runTime.deltaTValue(),2);
}

const bool eventTimeTracking =
    runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

// ************************************************************************* //
