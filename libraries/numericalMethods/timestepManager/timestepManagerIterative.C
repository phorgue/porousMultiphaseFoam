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

\*---------------------------------------------------------------------------*/

#include <label.H>
#include "timestepManagerIterative.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //
Foam::timestepManagerIterative::timestepManagerIterative(
    const Time& runTime,
    const dictionary& solutionDict,
    const word& algoName,
    const bool& steady
)
    :
    runTime_(runTime),
    iter_(0),
    iterIncrease_(-1),
    maxIter_(solutionDict.subOrEmptyDict(algoName).getOrDefault<label>("maxIter",10)),
    tolerance_(solutionDict.subOrEmptyDict(algoName).getOrDefault<scalar>("tolerance",1e+9)),
    nIterIncreasePresent_(runTime_.controlDict().found(("nIter"+algoName))),
    nIterIncrease_(runTime_.controlDict().getOrDefault<label>("nIter"+algoName, 0)),
    dTFactDecrease_(runTime_.controlDict().getOrDefault<scalar>("dTFactDecrease", 0.8)),
    dTFactIncrease_(runTime_.controlDict().getOrDefault<scalar>("dTFactIncrease", 1.25))
{
    if (steady)
    {
        maxIter_ = 1;
        if (algoName == "Newton") tolerance_ = 1e+9;
    }
    Info << nl << algoName << " loop control" << nl << "{"
    << nl << "    tolerance = " << tolerance_
    << nl << "    maximum number of iteration = " << maxIter_ << endl;

    if (nIterIncreasePresent_) Info << "    Number of iteration expected (for timestep increase) = " << nIterIncrease_ << endl;

    Info << "}" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timestepManagerIterative::~timestepManagerIterative()
{}

// * * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * //

scalar timestepManagerIterative::computeTimestep()
{
    scalar dt = runTime_.deltaTValue();
    if (iter_ > maxIter_) {
        dt = dTFactDecrease_ * runTime_.deltaTValue();
    }
    else {
        if (nIterIncreasePresent_) {
            if (iter_ < nIterIncrease_) {
                iterIncrease_++;
                if (iterIncrease_ == 5) {
                    dt = dTFactIncrease_ * runTime_.deltaTValue();
                    iterIncrease_ = 0;
                }
            }
            else {
                iterIncrease_ = 0;
            }
        }
        else {
            dt = GREAT;
        }

    }
    return dt;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
