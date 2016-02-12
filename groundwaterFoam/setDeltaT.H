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
    Timestep for Richards solver (Picard loop) 

\*---------------------------------------------------------------------------*/

if(adjustTimeStep)
{

    scalar deltaTFact=1;
    if (iterPicard > maxIterPicard)
    {
        deltaTFact = dTFactDecrease;
        iterStability=0;
    }
    else if (iterPicard < minIterPicard)
    {
        iterStability++;
    }

    if (iterStability == maxIterStability)
    {
        deltaTFact = dTFactIncrease;      
        iterStability=0;
    }
        
    //- Computing timestep 
    runTime.setDeltaT
        (
            min(
                deltaTFact*runTime.deltaTValue(),
                maxDeltaT
            )
        );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;

}

// ************************************************************************* //
