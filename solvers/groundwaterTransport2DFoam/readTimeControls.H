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

//- global parameters
bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
const bool eventTimeTracking =
     runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

//- water time step control
scalar dhmax =
    runTime.controlDict().lookupOrDefault<scalar>("dhmax",1e-6);

scalar dhdTmax = dhmax/runTime.deltaTValue();

//- tracer timestep control
scalar dCresidual =
    runTime.controlDict().lookupOrDefault<scalar>("dCresidual",1e-6);

scalar dCmax =
    runTime.controlDict().lookupOrDefault<scalar>("dCmax",0.01);

scalarList dCdTmax(composition.Y().size());
forAll(composition.Y(), speciesi)
{
    const auto& C = composition.Y(speciesi);
    dCdTmax[speciesi] = max(dCmax*gMax(C),dCresidual)/runTime.deltaTValue();
}

// ************************************************************************* //
