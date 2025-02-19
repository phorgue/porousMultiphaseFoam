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

#include "dualDynamicMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "dynamicRefineFvMesh.H"
#include "processorPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "linear.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dualDynamicMesh, 0);

addToRunTimeSelectionTable
(
        multiMesh,
        dualDynamicMesh,
        dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dualDynamicMesh::dualDynamicMesh
(
    Time& runTime,
    dynamicFvMesh& mesh
)
    :
    multiMesh(runTime, mesh),
    fineMeshPtr_(nullptr),
    scalarFields_(),
    vectorFields_(),
    phiFields_()
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
                if (not(isA<processorPolyPatch>(patches[patchi])) && not(isA<symmetryPolyPatch>(patches[patchi]))){
                    boundary_protected_cells.append(fineMesh.boundary()[patchi].faceCells());
                }
            }
            refCast<dynamicRefineFvMesh>(fineMesh).protectedCell() = bitSet(boundary_protected_cells);
        }
    }
}

// * * * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField>
void Foam::dualDynamicMesh::mapFieldCoarseToFine(
    Foam::GeometricField<Type, PatchField, Foam::volMesh>& field1,
    Foam::GeometricField<Type, PatchField, Foam::volMesh>& field2
)
{
    dynamicFvMesh& fineMesh = fineMeshPtr_.ref();
    const refinementHistory& ref_hist = refCast<dynamicRefineFvMesh>(fineMesh).meshCutter().history();
    const labelIOList& cell_level =  refCast<dynamicRefineFvMesh>(fineMesh).meshCutter().cellLevel();
    label list_size = ref_hist.splitCells().size();

    List<Type> savedValues(list_size, Foam::Zero);
    List<Type> savedValuesOld(list_size, Foam::Zero);

    for(label celli=0;celli<coarseMesh_.nCells();celli++)
    {
        field2[celli] = field1[celli];
        field2.oldTime()[celli] = field1.oldTime()[celli];

        // propagation of value
        if (cell_level[celli] > 0)
        {
            label cell0 = ref_hist.splitCells()[ref_hist.visibleCells()[celli]].parent_;
            for(label iter=1;iter<cell_level[celli];iter++)
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
            label cell0 = ref_hist.splitCells()[ref_hist.visibleCells()[celli]].parent_;
            for(label iter=1;iter<cell_level[celli];iter++)
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
        if (isA<processorPolyPatch>(fineMesh.boundaryMesh()[patchi]) ||
            isA<symmetryPolyPatch>(fineMesh.boundaryMesh()[patchi]) )
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
    field2.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::dualDynamicMesh::addField
(
    volScalarField& coarseField
)
{
    volScalarField* coarsePointer = &coarseField;
    scalarFields_.first().append(coarsePointer);
    label current_id = scalarFields_.first().size()-1;
    Info << nl << "create dual scalar field " << current_id << "  for " << coarseField.name() << "...";
    scalarFields_.second().append(new volScalarField(
            IOobject
            (
                coarseField.name()+"_dual",
                coarseField.time().timeName(),
                fineMeshPtr_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fineMeshPtr_,
            0,
            coarseField.dimensions()
        )
    );
    mapFieldCoarseToFine(coarseField, *scalarFields_.second()[current_id]);
    Info << "ok" << endl;
    return *scalarFields_.second()[current_id];
}

Foam::volVectorField& Foam::dualDynamicMesh::addField
    (
        volVectorField& coarseField
    )
{
    volVectorField* coarsePointer = &coarseField;
    vectorFields_.first().append(coarsePointer);
    label current_id = vectorFields_.first().size()-1;
    Info << nl << "create dual vector field " << current_id << "  for " << coarseField.name() << "...";
    vectorFields_.second().append(new volVectorField(
            IOobject
            (
                coarseField.name()+"_dual",
                coarseField.time().timeName(),
                fineMeshPtr_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fineMeshPtr_,
            Vector<scalar>(0,0,0),
            coarseField.dimensions()
        )
    );
    mapFieldCoarseToFine(coarseField, *vectorFields_.second()[current_id]);
    Info << "ok" << endl;
    return *vectorFields_.second()[current_id];
}

Foam::surfaceScalarField& Foam::dualDynamicMesh::addField
        (
                surfaceScalarField& coarsePhi
        )
{
    Info << nl << "Set dual flux field...";
    phiFields_.first() = &coarsePhi;
    phiFields_.second() = new surfaceScalarField
    (
        IOobject
            (
                coarsePhi.name()+"_dual",
                coarsePhi.time().timeName(),
                fineMeshPtr_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
        fineMeshPtr_,
        0,
        coarsePhi.dimensions(),
        coarsePhi.boundaryField().types()
    );
    Info << "ok" << endl;
    return *phiFields_.second();
}


bool Foam::dualDynamicMesh::update()
{
    fineMeshPtr_->update();
    for(label i=0;i<scalarFields_.first().size();i++) {
        mapFieldCoarseToFine(*scalarFields_.first()[i], *scalarFields_.second()[i]);
    }
    for(label i=0;i<vectorFields_.first().size();i++) {
        mapFieldCoarseToFine(*vectorFields_.first()[i], *vectorFields_.second()[i]);
    }
    volVectorField& vField = *vectorFields_.second()[0];
    *phiFields_.second() = linearInterpolate(vField) & vField.mesh().Sf();

    if (fineMeshPtr_->changing()) return true;
    else return false;
}

// ************************************************************************* //
