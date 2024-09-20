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
    volScalarField& coarseField,
    volScalarField& fineField
)
{
    volScalarField* coarsePointer = &coarseField;
    volScalarField* finePointer = &fineField;
    Tuple2<volScalarField*, volScalarField*> tfields(coarsePointer, finePointer);
    scalarFields_.append(tfields);
}

void Foam::dualMesh::addDualFields
(
    volVectorField& coarseField,
    volVectorField& fineField
)
{
    volVectorField* coarsePointer = &coarseField;
    volVectorField* finePointer = &fineField;
    Tuple2<volVectorField*, volVectorField*> tfields(coarsePointer, finePointer);
    vectorFields_.append(tfields);
}

void Foam::dualMesh::update()
{
    dynamicFvMesh& fineMesh = fineMeshPtr_.ref();
    const refinementHistory& ref_hist = refCast<dynamicRefineFvMesh>(fineMesh).meshCutter().history();
    const labelIOList& cell_level =  refCast<dynamicRefineFvMesh>(fineMesh).meshCutter().cellLevel();
    label list_size = ref_hist.splitCells().size();

    label cell0 = 0;
    forAll(scalarFields_, fieldsi)
    {
        volScalarField& field1 = *scalarFields_[fieldsi].first();
        volScalarField& field2 = *scalarFields_[fieldsi].second();
        scalarList savedValues(list_size, 0);
        scalarList savedValuesOld(list_size, 0);


        for(label celli=0;celli<coarseMesh_.nCells();celli++)
        {
            field2[celli] = field1[celli];
            field2.oldTime()[celli] = field1.oldTime()[celli];

            // propagation of value
            if (cell_level[celli] > 0)
            {
                for(label iter=0;iter<cell_level[celli];iter++)
                {
                    cell0 = ref_hist.splitCells()[cell0].parent_;
                }
                cell0 /= 9;
                savedValues[cell0] = field1[celli];
                savedValuesOld[cell0] = field1.oldTime()[celli];
            }
        }

        if (coarseMesh_.nCells() != fineMesh.nCells())
        {
            // labels of all new cells created by the remeshing
            for(label celli=coarseMesh_.nCells();celli<fineMesh.nCells();celli++)
            {
                for(label iter=0;iter<cell_level[celli];iter++)
                {
                    cell0 = ref_hist.splitCells()[cell0].parent_;
                }
                cell0 /= 9;
                field2[celli] = savedValues[cell0];
                field2.oldTime()[celli] = savedValuesOld[cell0];
            }
        }

        forAll(fineMesh.boundary(), patchi)
        {
            if (isA<processorPolyPatch>(fineMesh.boundaryMesh()[patchi]))
            {
                field2.boundaryFieldRef()[patchi].initEvaluate(Pstream::commsTypes::nonBlocking);
                field2.oldTime().boundaryFieldRef()[patchi].initEvaluate(Pstream::commsTypes::nonBlocking);
            }
            else
            {
                forAll(field2.boundaryFieldRef()[patchi], facei)
                {
                    field2.boundaryFieldRef()[patchi][facei] = field1.boundaryField()[patchi][facei];
                    field2.oldTime().boundaryFieldRef()[patchi][facei] = field1.oldTime().boundaryField()[patchi][facei];
                }
            }
        }

    }

}

// ************************************************************************* //
