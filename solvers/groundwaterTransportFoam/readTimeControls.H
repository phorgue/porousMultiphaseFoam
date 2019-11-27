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
    Read the control parameters used by setDeltaT for groundwaterTransportFoam
    solver

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
const scalar dTFactIncrease = runTime.controlDict().lookupOrDefault<scalar>("dTFactIncrease",0);
const scalar dTFactDecrease = runTime.controlDict().lookupOrDefault<scalar>("dTFactDecrease",0);

//- if dtFromtheta used
const scalar dthetamax = runTime.controlDict().lookupOrDefault("dthetamax", 1e-2);
scalar dthetadTmax = dthetamax/runTime.deltaTValue();
scalar dtheta = 0;
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
            << nl << "    dthetamax   = " << dthetamax
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

//- tracer timestep control
scalar truncationError =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError",0.01);
word timeScheme = mesh.ddtScheme("ddt(" + composition.Y(0).name() + ')');
EulerD3dt3Scheme d3dt3Operator(mesh,runTime.deltaTValue());
Foam::fv::EulerD2dt2Scheme<scalar> d2dt2Operator(mesh);
scalarList dCdTmax(composition.Y().size());
scalarList dC2dT2max(composition.Y().size());
scalarList dC3dT3max(composition.Y().size());
scalarList Cmax(composition.Y().size());
//- Derivative initialisation to keep 1st user-defined time step
if (timeScheme == "backward")
{
    forAll(composition.Y(), speciesi)
    {
        const auto& C = composition.Y(speciesi);
        Cmax[speciesi] = gMax(C);
        dC3dT3max[speciesi] = 3*truncationError*(Cmax[speciesi]+VSMALL)/Foam::pow(runTime.deltaTValue(),3);
    }
}
else if  (timeScheme == "CrankNicolson")
{
    forAll(composition.Y(), speciesi)
    {
        const auto& C = composition.Y(speciesi);
        Cmax[speciesi] = gMax(C);
        dC3dT3max[speciesi] = 12*truncationError*(Cmax[speciesi]+VSMALL)/Foam::pow(runTime.deltaTValue(),3);
    }
}
else if (timeScheme == "Euler")
{
    forAll(composition.Y(), speciesi)
    {
        const auto& C = composition.Y(speciesi);
        Cmax[speciesi] = gMax(C);
        dC2dT2max[speciesi] = 2*truncationError*(Cmax[speciesi]+VSMALL)/Foam::pow(runTime.deltaTValue(),2);
    }
}

// ************************************************************************* //
