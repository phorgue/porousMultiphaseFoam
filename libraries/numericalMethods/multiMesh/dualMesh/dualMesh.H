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
    Foam::dualMesh

Description
    A class for handling fields transfert between dual meshes.

SourceFiles
    dualMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dualMesh_H
#define dualMesh_H

#include "multiMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class porousMediumModel Declaration
\*---------------------------------------------------------------------------*/

class dualMesh
:
    public multiMesh
{
protected:

    // Protected data
    //- Reference to secondary mesh
    autoPtr<dynamicFvMesh> fineMeshPtr_;

    //- List of tuple_field
    List<Tuple2<const volScalarField*, volScalarField*>> scalarFields_;
    List<Tuple2<const volVectorField*, volVectorField*>> vectorFields_;

    //- Disallow copy construct
    dualMesh(const dualMesh&);

    //- Disallow default bitwise assignment
    void operator=(const dualMesh&);

public:

    //- Runtime type information
    TypeName("dualMesh");

    //- Construct from components
    dualMesh(dynamicFvMesh& mesh);

    //- Destructor
    virtual ~dualMesh() {}

    //- Return mesh
    dynamicFvMesh& fineMesh() {return fineMeshPtr_.ref();}

    //- Increment list of dual-mesh fields
    void addDualFields(const volScalarField&, volScalarField&);
    void addDualFields(const volVectorField&, volVectorField&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
