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
    Foam::dualDynamicMesh

Description
    A class for handling fields transfert between dual meshes.

SourceFiles
    dualDynamicMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dualDynamicMesh_H
#define dualDynamicMesh_H

#include "multiMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class porousMediumModel Declaration
\*---------------------------------------------------------------------------*/

class dualDynamicMesh
:
    public multiMesh
{
protected:

    // Protected data
    //- Reference to secondary mesh
    autoPtr<dynamicFvMesh> fineMeshPtr_;

    //- List of tuple_field
    Tuple2<PtrList<volScalarField>, PtrList<volScalarField>> scalarFields_;
    Tuple2<PtrList<volVectorField>, PtrList<volVectorField>> vectorFields_;
    PtrList<surfaceScalarField> phiFields_;

    //- Function to map fields between dual meshes
    template<class Type, template<class> class PatchField>
        void mapFieldCoarseToFine(
            GeometricField<Type, PatchField, volMesh>& field1,
            GeometricField<Type, PatchField, volMesh>& field2
            );

    //- Disallow copy construct
    dualDynamicMesh(const dualDynamicMesh&);

    //- Disallow default bitwise assignment
    void operator=(const dualDynamicMesh&);

public:

    //- Runtime type information
    TypeName("dualDynamicMesh");

    //- Construct from components
    dualDynamicMesh(dynamicFvMesh& mesh);

    //- Destructor
    virtual ~dualDynamicMesh() {}

    //- Return mesh
    dynamicFvMesh& fineMesh() {return fineMeshPtr_.ref();}

    //- Increment list of dual-mesh fields
    volScalarField& addField(volScalarField&);
    volVectorField& addField(volVectorField&);
    surfaceScalarField& addField(surfaceScalarField&);

    bool dynamic() { return true;}
    void update();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
