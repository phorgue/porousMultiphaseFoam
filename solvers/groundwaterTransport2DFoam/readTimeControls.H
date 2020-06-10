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

//- tracer and hwater timestep control
scalar truncationError_global =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError",0.01);

scalar truncationError_potential =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError_potential",truncationError_global);

scalar truncationError_C =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError_C",truncationError_global);

word timeScheme = mesh.ddtScheme("ddt(" + composition.Y(0).name() + ')');
EulerD3dt3Scheme d3dt3Operator(mesh,runTime.deltaTValue());
Foam::fv::EulerD2dt2Scheme<scalar> d2dt2Operator(mesh);
scalar dpotentialdTmax = 0;
scalar dpotential2dT2max = 0;
scalar dpotential3dT3max = 0;
scalar hwatermax = 0;
scalarList dCdTmax(composition.Y().size());
scalarList dC2dT2max(composition.Y().size());
scalarList dC3dT3max(composition.Y().size());
scalarList Cmax(composition.Y().size());

if (adjustTimeStep)
{
    //- Derivative initialisation to keep 1st user-defined time step
    if (timeScheme == "backward")
    {
        forAll(composition.Y(), speciesi)
        {
            const auto& C = composition.Y(speciesi);
            Cmax[speciesi] = gMax(C);
            dC3dT3max[speciesi] = 3*truncationError_C*(Cmax[speciesi]+VSMALL)/Foam::pow(runTime.deltaTValue(),3);
        }
        hwatermax = gMax(hwater);
        dpotential3dT3max = 3*truncationError_potential*(hwatermax+VSMALL)/Foam::pow(runTime.deltaTValue(),3);
    }
    else if  (timeScheme == "CrankNicolson")
    {
        forAll(composition.Y(), speciesi)
        {
            const auto& C = composition.Y(speciesi);
            Cmax[speciesi] = gMax(C);
            dC3dT3max[speciesi] = 12*truncationError_C*(Cmax[speciesi]+VSMALL)/Foam::pow(runTime.deltaTValue(),3);
        }
        hwatermax = gMax(hwater);
        dpotential3dT3max = 12*truncationError_potential*(hwatermax+VSMALL)/Foam::pow(runTime.deltaTValue(),3);
    }
    else if (timeScheme == "Euler")
    {
        forAll(composition.Y(), speciesi)
        {
            const auto& C = composition.Y(speciesi);
            Cmax[speciesi] = gMax(C);
            dC2dT2max[speciesi] = 2*truncationError_C*(Cmax[speciesi]+VSMALL)/Foam::pow(runTime.deltaTValue(),2);
        }
        hwatermax = gMax(hwater);
        dpotential2dT2max = 2*truncationError_potential*(hwatermax+VSMALL)/Foam::pow(runTime.deltaTValue(),2);
    }
    Info << "Timestepping for potential equation is based on time-scheme truncation error with :"
        << nl << "{"
        << nl << "    truncationError = " << truncationError_potential
        << nl << "}"
        << endl;
    Info << "Timestepping for C equation is based on time-scheme truncation error with :"
        << nl << "{"
        << nl << "    truncationError = " << truncationError_C
        << nl << "}"
        << endl;
}

const bool eventTimeTracking =
     runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

// ************************************************************************* //
