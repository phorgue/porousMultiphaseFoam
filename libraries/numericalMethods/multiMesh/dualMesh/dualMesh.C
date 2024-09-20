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

#include "dualMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "dynamicRefineFvMesh.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dualMesh, 0);

addToRunTimeSelectionTable
(
    multiMesh,
    dualMesh,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dualMesh::dualMesh
(
    dynamicFvMesh& mesh
)
    :
    multiMesh(mesh),
    fineMeshPtr_(nullptr),
    scalarFields_(0),
    vectorFields_(0)
{
    //- Read the dual mesh
    fineMeshPtr_ = dynamicFvMesh::New
    (
        IOobject
        (
            "dualMesh",
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ
        )
    );

    //- Protect patch dual-mesh cells from refinement
    dynamicFvMesh& fineMesh = fineMeshPtr_.ref();
    if (fineMesh.dynamic()) {
        if (isA<dynamicRefineFvMesh>(fineMesh)) {
            DynamicList<label> boundary_protected_cells(0);
            const polyBoundaryMesh& patches = fineMesh.boundaryMesh();
            forAll(patches, patchi) {
                if (not(isA<processorPolyPatch>(patches[patchi])) &&
                    not(isA<symmetryPlanePolyPatch>(patches[patchi])) ){
                    boundary_protected_cells.append(fineMesh.boundary()[patchi].faceCells());
                }
            }
            refCast<dynamicRefineFvMesh>(fineMesh).protectedCell() = bitSet(boundary_protected_cells);
        }
    }
}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

void Foam::dualMesh::addDualFields
(
    const volScalarField& coarseField,
    volScalarField& fineField
)
{
    const volScalarField* coarsePointer = &coarseField;
    volScalarField* finePointer = &fineField;
    Tuple2<const volScalarField*, volScalarField*> tfields(coarsePointer, finePointer);
    scalarFields_.append(tfields);
}

void Foam::dualMesh::addDualFields
(
    const volVectorField& coarseField,
    volVectorField& fineField
)
{
    const volVectorField* coarsePointer = &coarseField;
    volVectorField* finePointer = &fineField;
    Tuple2<const volVectorField*, volVectorField*> tfields(coarsePointer, finePointer);
    vectorFields_.append(tfields);
}

// ************************************************************************* //
