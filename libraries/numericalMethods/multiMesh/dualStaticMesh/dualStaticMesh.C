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

#include "dualStaticMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "linear.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dualStaticMesh, 0);

addToRunTimeSelectionTable
(
    multiMesh,
    dualStaticMesh,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dualStaticMesh::dualStaticMesh
(
    Time& runTime,
    dynamicFvMesh& mesh
)
    :
    multiMesh(runTime, mesh),
    fineMeshPtr_(nullptr),
    runTime_(runTime),
    cellMapping_(0),
    boundaryMapping_(0),
    refined_(false),
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
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
}

// * * * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

void Foam::dualStaticMesh::initialRefinement()
{
    //- Refining mesh
    Info << nl << "**** Initialization of static dual mesh *****" << endl;
    dynamicFvMesh& fineMesh = fineMeshPtr_();
    label timeIndex = runTime_.timeIndex();
    runTime_.setTime(runTime_.timeOutputValue(), timeIndex+1);
    fineMesh.movePoints(fineMesh.points());
    label ncells = 0;
    while (fineMesh.nCells() != ncells) {
        ncells = fineMesh.nCells();
        fineMesh.update();
    }
    runTime_.setTime(runTime_.timeOutputValue(), timeIndex);

    //- Contruct cell mapping
    cellMapping_.resize(fineMesh.nCells(), -1);
    Info << "Construct cell mapping between coarse and fine mesh...";
    forAll(fineMesh.cells(), tgtCelli) {
        const point& tgtPos = fineMesh.cellCentres()[tgtCelli];
        label foundCell = coarseMesh_.findCell(tgtPos);
        if (foundCell != -1){
            cellMapping_[tgtCelli] = foundCell;
        }
        else
        {
            FatalErrorIn("dualStaticMesh.C") << "dual mesh addressing error, decomposition of"
            "coarse and fine mesh should be the same" << abort(FatalError);
        }
    }
    Info << "ok" << endl;

    //-Construct boundary mapping
    const polyBoundaryMesh& bCoarseMesh = coarseMesh_.boundaryMesh();
    const polyBoundaryMesh& bFineMesh = fineMesh.boundaryMesh();
    const surfaceVectorField& cSf = coarseMesh_.Sf();
    const surfaceVectorField& fSf = fineMesh.Sf();

    boundaryMapping_.resize(bFineMesh.size());
    forAll(fSf.boundaryField(), patchi) {
        boundaryMapping_[patchi].resize(fSf.boundaryField()[patchi].size());
        forAll(fSf.boundaryField()[patchi], tgtFacei) {
            label fCell = fineMesh.boundary().faceCells()[patchi][tgtFacei];
            label cCell = cellMapping_[fCell];
            forAll(coarseMesh_.cells()[cCell], facei) {
                label b_face = coarseMesh_.cells()[cCell][facei] - bCoarseMesh.patchStarts()[patchi];
                if (b_face > -1 && b_face < bCoarseMesh.patchSizes()[patchi]) {
                    if ((fSf.boundaryField()[patchi][tgtFacei] & cSf.boundaryField()[patchi][b_face]) > 0 ) {
                        boundaryMapping_[patchi][tgtFacei] = b_face;
                    }
                }
            }
        }

    }

    for(label i=0;i<scalarFields_.first().size();i++) {
        mapFieldCoarseToFine(*scalarFields_.first()[i], *scalarFields_.second()[i]);
    }
    for(label i=0;i<vectorFields_.first().size();i++) {
        mapFieldCoarseToFine(*vectorFields_.first()[i], *vectorFields_.second()[i]);
    }
    Info << "*********************************************" << endl << endl;
}

template<class Type, template<class> class PatchField>
void Foam::dualStaticMesh::mapFieldCoarseToFine(
    Foam::GeometricField<Type, PatchField, Foam::volMesh>& field1,
    Foam::GeometricField<Type, PatchField, Foam::volMesh>& field2
)
{
    forAll(field2, celli) {
        field2[celli] = field1[cellMapping_[celli]];
        field2.oldTime()[celli] = field1.oldTime()[cellMapping_[celli]];
    }
    forAll(field2.boundaryField(), patchi) {
        forAll(field2.boundaryField()[patchi], facei) {
            field2.boundaryFieldRef()[patchi][facei] =
                    field1.boundaryField()[patchi][boundaryMapping_[patchi][facei]];
            field2.oldTime().boundaryFieldRef()[patchi][facei] =
                    field1.oldTime().boundaryField()[patchi][boundaryMapping_[patchi][facei]];
        }
    }
}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::dualStaticMesh::addField
(
    volScalarField& coarseField
)
{
    volScalarField* coarsePointer = &coarseField;
    scalarFields_.first().append(coarsePointer);
    label current_id = scalarFields_.first().size()-1;
    Info << nl << "Create dual scalar field " << current_id << " for " << coarseField.name() << "...";
    scalarFields_.second().append(new volScalarField(
            IOobject
            (
                coarseField.name()+"_dual",
                coarseField.time().timeName(),
                fineMeshPtr_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fineMeshPtr_,
            0,
            coarseField.dimensions(),
            coarseField.boundaryField().types()
        )
    );
    scalarFields_.second()[current_id]->primitiveFieldRef() = coarseField.primitiveField();
    forAll(coarseField.boundaryField(), patchi)
    {
        forAll(coarseField.boundaryField()[patchi], facei)
        {
            scalarFields_.second()[current_id]->boundaryFieldRef()[patchi][facei] =
                coarseField.boundaryField()[patchi][facei];
        }
    }
    Info << "ok" << endl;
    scalarFields_.second()[current_id]->write();
    return *scalarFields_.second()[current_id];
}

Foam::volVectorField& Foam::dualStaticMesh::addField
    (
        volVectorField& coarseField
    )
{
    volVectorField* coarsePointer = &coarseField;
    vectorFields_.first().append(coarsePointer);
    label current_id = vectorFields_.first().size()-1;
    Info << nl << "create dual vector field " << current_id << " for " << coarseField.name() << "...";
    vectorFields_.second().append(new volVectorField(
            IOobject
            (
                coarseField.name()+"_dual",
                coarseField.time().timeName(),
                fineMeshPtr_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fineMeshPtr_,
            Vector<scalar>(0,0,0),
            coarseField.dimensions(),
            coarseField.boundaryField().types()
        )
    );
    vectorFields_.second()[current_id]->primitiveFieldRef() = coarseField.primitiveField();
    vectorFields_.second()[current_id]->correctBoundaryConditions();
    Info << "ok" << endl;
    vectorFields_.second()[current_id]->write();
    return *vectorFields_.second()[current_id];
}

Foam::surfaceScalarField& Foam::dualStaticMesh::addField
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
            IOobject::NO_READ,
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

bool Foam::dualStaticMesh::update()
{
    bool meshChanged = false;
    if (!refined_)
    {
        initialRefinement();
        refined_ = true;
        meshChanged = true;
    }
    for(label i=0;i<scalarFields_.first().size();i++) {
        mapFieldCoarseToFine(*scalarFields_.first()[i], *scalarFields_.second()[i]);
    }
    for(label i=0;i<vectorFields_.first().size();i++) {
        mapFieldCoarseToFine(*vectorFields_.first()[i], *vectorFields_.second()[i]);
    }

    volVectorField& vField = *vectorFields_.second()[0];
    *phiFields_.second() = linearInterpolate(vField) & vField.mesh().Sf();

    return meshChanged;
}

// ************************************************************************* //
