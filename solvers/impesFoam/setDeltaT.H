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

Description
    Timestep for IMPES solver 

\*---------------------------------------------------------------------------*/

//- Timestep Saturation equation
scalarField dSdt= mag(-fvc::div(phib)+sourceTerm)().internalField()/eps.internalField();
scalar dSdTmax = gMax(dSdt);
scalar dtForS = dSmax/(dSdTmax+SMALL);

//- Choose between Coats and Courant number as driving number
scalar CFLUse = -1;
scalar maxDeltaTFact = -1;
switch(static_cast<int>(CFLint))
{
    case 0:
        #include "CoatsNo.H"
        break;
    case 1:
        #include "ToddNo.H"
        break;
    case 2:
        #include "CourantNo.H"
        break;
};

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
