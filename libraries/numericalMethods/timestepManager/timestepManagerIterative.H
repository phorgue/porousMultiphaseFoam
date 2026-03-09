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

Class
    Foam::fv::timestepManagerIterative

Description
    Timestep manager compute timestep based on number of linear algorithm's
    iterations.

SourceFiles
    timestepManagerIterative.C

\*---------------------------------------------------------------------------*/

#ifndef timestepManagerIterative_H
#define timestepManagerIterative_H

#include "nonLinearEqn.H"
#include "fvMesh.H"
#include "EulerD2dt2Scheme.H"
#include "EulerD3dt3Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Class timestepManagerIterative Declaration
\*---------------------------------------------------------------------------*/

class timestepManagerIterative
{
    
private:

    //- Reference to time
    const Time& runTime_;

    const word algoName_;

    //- Solver iterator
    label iter_;
    label iterIncrease_;

    //- Parameters
    label maxIter_;
    scalar tolerance_;
    const bool nIterIncreasePresent_;
    const label nIterIncrease_;
    const scalar dTFactDecrease_;
    const scalar dTFactIncrease_;

    Tuple2<scalar, scalar> residualDelta_;

public:

    // Constructors

    //- Construct from field
    timestepManagerIterative(
        const Time&,
        const dictionary&,
        const word&,
        const bool&
    );

    virtual ~timestepManagerIterative();

    //- Access function
    void reset(const label value = 0) {iter_ = value;}
    const label& iter() const {return iter_;}
    const label& maxIter() const {return maxIter_;}
    const scalar& tolerance() const {return tolerance_;}
    const Tuple2<scalar, scalar>& residualDelta() const {return residualDelta_;}

    //- Member function
    void operator++(int) {iter_++;};
    scalar computeTimestep();
    bool solveEquation(nonLinearEqn&, bool, label, scalar initResidual = 1.0001);
    

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
