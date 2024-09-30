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
    Foam::multiMesh

Description
    An abstract base class for simple or dual meshes.

SourceFiles
    multiMesh.C
    multiMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef multiMesh_H
#define multiMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class multiMesh Declaration
\*---------------------------------------------------------------------------*/

class multiMesh
{
protected:

    //- Reference to the main/coarse mesh
    dynamicFvMesh& coarseMesh_;

public:

    //- Runtime type information
    TypeName("multiMesh");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        multiMesh,
        dictionary,
        (
        dynamicFvMesh& mesh
        ),
    (mesh)
    );

    // Selectors

    //- Return a reference to the selected multi mash model
    static autoPtr<multiMesh> New
    (
        dynamicFvMesh& mesh,
        const dictionary& transportProperties
    );

    // Constructors

    //- Construct from components
    multiMesh(dynamicFvMesh& mesh);

    //- Destructor
    virtual ~multiMesh()
        {}

    //- Return mesh
    dynamicFvMesh& coarseMesh() {return coarseMesh_;}

    //- Return transport properties
    virtual dynamicFvMesh& fineMesh() = 0;

    //- Add base and dual mesh fields
    virtual volScalarField& addField(volScalarField&) = 0;
    virtual volVectorField& addField(volVectorField&) = 0;
    virtual surfaceScalarField& addField(surfaceScalarField&) = 0;
    virtual void update() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
