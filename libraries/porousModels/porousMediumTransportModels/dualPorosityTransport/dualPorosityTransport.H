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

Namespace
    Foam::porousMediumModels::dualPorosityTransport

Class
    Foam::dualPorosityTransport

Description
    Dual porosity model (one kr / one pc law)

SourceFiles
    dualPorosityTransport.C

\*---------------------------------------------------------------------------*/

#ifndef dualPorosityTransport_H
#define dualPorosityTransport_H

#include "volFields.H"
#include "surfaceFields.H"
#include "porousMediumModelTransport.H"
#include "incompressiblePhase.H"
#include "mul

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumTransportModels
{

/*---------------------------------------------------------------------------* \
             Class dualPorosityTransport Declaration
\*---------------------------------------------------------------------------*/

class dualPorosityTransport
:
    public porousMediumTransportModel
{
    const dictionary dualPorosityTransportCoeffs_;
    const fvMesh& mesh_;

    //- Disallow copy construct
    dualPorosityTransport(const dualPorosityTransport&);

    //- Disallow default bitwise assignment
    void operator=(const dualPorosityTransport&);

public:

    //- Runtime type information
    TypeName("dualPorosityTransport");

    // Constructors

    //- Construct from components
    dualPorosityTransport
    (
        const word Sname,
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const porousMediumModel& pmModel
)
    );

    //- Destructor
    ~dualPorosityTransport()
        {}

    //- Correct
    void correct();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousMediumTransportModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
