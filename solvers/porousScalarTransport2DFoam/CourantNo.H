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
    
Global
     CourantNo
    
Description
     Calculates and outputs the mean and maximum Courant Numbers.
    
\*---------------------------------------------------------------------------*/
    
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
    
{
    scalarField sumPhi
        (
            fvc::surfaceSum(mag(phihwater))().primitiveField()
        );
    
    CoNum = 0.5*gMax(sumPhi/(eps.field()*hwater.field()*mesh.V().field()))*runTime.deltaTValue();
    
    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(eps.field()*hwater.field()*mesh.V().field()))*runTime.deltaTValue();
}
    
Info<<  "Courant Number mean: " << meanCoNum << " max: " << CoNum << endl;
    
// ************************************************************************* //
