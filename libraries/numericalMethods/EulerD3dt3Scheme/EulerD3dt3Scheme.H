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
    Foam::fv::EulerD3dt3Scheme

Description
    First-order Euler explicit d3dt3 using the current and two previous
    time-step values.

SourceFiles
    EulerD3dt3Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef EulerD3dt3Scheme_H
#define EulerD3dt3Scheme_H

#include "fvMesh.H"
#include "tmp.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Class EulerD3dt3Scheme Declaration
\*---------------------------------------------------------------------------*/

class EulerD3dt3Scheme
{
    
protected:
   
    const fvMesh& mesh_; 

    //- store old-old-time step
    scalar dt00_;
    
public:

    // Constructors

        //- Construct from mesh
    EulerD3dt3Scheme(const fvMesh& mesh, scalar dt)
        :
            mesh_(mesh),
            dt00_(dt)
        {}

    virtual ~EulerD3dt3Scheme();

    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        tmp<volScalarField> d3dt3
        (
            const volScalarField&
        );

    // Member Operators
       //- Store the old-old time step
        void storeDeltaT00(scalar dt){ dt00_ = dt;}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
