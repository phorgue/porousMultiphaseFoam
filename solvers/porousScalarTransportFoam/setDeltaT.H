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
 
 Global
     setDeltaT
 
 Description
     Reset the timestep to depending of C variation between two timesteps.
     Reduction of time-step is immediate, but increase is limited to 20% to
     avoid unstable oscillations.
 
 \*---------------------------------------------------------------------------*/
   
 if (adjustTimeStep)
 {
     scalar dCmax = max(dCmax_rel*gMax(C),dCmax_abs);

     scalar dtForC =  dCmax/(dCdTmax+SMALL);
     scalar newDeltaT = min(dtForC, 1.2*runTime.deltaTValue());
 
     runTime.setDeltaT
     (
         min
         (
             newDeltaT,
             maxDeltaT
         )
     );
 
     Info<< nl << "deltaT = " <<  runTime.deltaTValue() << " dCmax = " << dCmax << endl;
 }
 
 // ************************************************************************* //
