/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Timestep for IMPES solver 

\*---------------------------------------------------------------------------*/

//- Timestep Saturation equation
scalar dSdTmax = gMax(mag(-fvc::div(phib)+(SrcInj*Winj-SrcExt*Wext)*activateWellbores)().internalField()/eps.internalField());
scalar dtForS = dSmax/(dSdTmax+SMALL);

//- Choose between Coats and Courant number as driving number
CFL =  runTime.controlDict().lookupOrDefault<bool>("CFL", true);
scalar CFLUse = CFL ? gMax(CFLCoats) : CoNum;

//- Caracteristic number
scalar maxDeltaTFact = maxCo/(CFLUse + SMALL);
scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

//- Computing timestep 
runTime.setDeltaT
(
    min(
        dtForS,
        min(
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    )
);

Info<< "deltaT = " <<  runTime.deltaTValue() << endl;

// ************************************************************************* //
