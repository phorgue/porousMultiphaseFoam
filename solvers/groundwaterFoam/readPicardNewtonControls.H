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
    readPicardNewtonControls

Description
    Read control for Picard and Newton's iterations

\*---------------------------------------------------------------------------*/

const dictionary& PicardDict = mesh.solutionDict().subOrEmptyDict("Picard");
const scalar tolerancePicard = PicardDict.lookupOrDefault<scalar>("tolerance", GREAT);
label maxIterPicard = PicardDict.lookupOrDefault<label>("maxIter",10);
const dictionary& NewtonDict = mesh.solutionDict().subOrEmptyDict("Newton");
const scalar toleranceNewton = NewtonDict.lookupOrDefault<scalar>("tolerance", GREAT);
const label maxIterNewton = NewtonDict.lookupOrDefault<label>("maxIter",10);

Info << nl << "Picard loop control" << nl << "{"
<< nl << "    tolerance = " << tolerancePicard
<< "    maximum number of iteration = " << maxIterPicard
<< nl << "}" << endl;
Info << nl << "Newton loop control" << nl << "{"
<< nl << "    tolerance = " << toleranceNewton
<< "    maximum number of iteration = " << maxIterNewton
<< nl << "}" << endl;

if (!PicardDict.found("tolerance") && !NewtonDict.found("tolerance"))
{
    FatalErrorIn("readPicardNewtonControls.H")
        << "Missing keyword tolerance in system/fvSolution.Picard or system/fvSolution.Newton"
            << abort(FatalError);
}

//- to switch between conservative mixed formulation (by default) and head pressure (non conservative)
bool massConservative = transportProperties.lookupOrDefault<bool>("massConservative",true);
