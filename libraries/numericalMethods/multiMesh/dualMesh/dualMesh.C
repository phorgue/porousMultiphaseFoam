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
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "linear.H"

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
    mapping_(nullptr),
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

    mapping_.reset(new meshToMesh(mesh, fineMeshPtr_, "mapNearest", "nearestFaceAMI"));
}

// * * * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField>
void Foam::dualMesh::mapFieldCoarseToFine(
    Foam::GeometricField<Type, PatchField, Foam::volMesh>& field1,
    Foam::GeometricField<Type, PatchField, Foam::volMesh>& field2
)
{
    field2 = mapping_.ref().mapSrcToTgt(field1);
}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::dualMesh::addField
(
    volScalarField& coarseField
)
{
    volScalarField* coarsePointer = &coarseField;
    scalarFields_.first().append(coarsePointer);
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
            coarseField.dimensions()
        )
    );
    mapFieldCoarseToFine(coarseField, scalarFields_.second().back());
    Info << nl << "create dual field for " << coarseField.name() << endl;
    scalarFields_.second().back().write();
    return scalarFields_.second().back();
}

Foam::volVectorField& Foam::dualMesh::addField
    (
        volVectorField& coarseField
    )
{
    volVectorField* coarsePointer = &coarseField;
    vectorFields_.first().append(coarsePointer);
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
            coarseField.dimensions()
        )
    );
    mapFieldCoarseToFine(coarseField, vectorFields_.second().back());
    Info << nl << "create dual field for " << coarseField.name() << endl;
    vectorFields_.second().back().write();
    return vectorFields_.second().back();
}

Foam::surfaceScalarField& Foam::dualMesh::addField
        (
                surfaceScalarField& coarseField
        )
{
    volVectorField& vField = vectorFields_.second().back();
    phiFields_.append(new surfaceScalarField
    (
        IOobject
            (
                coarseField.name()+"_dual",
                vField.time().timeName(),
                vField.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
        linearInterpolate(vField) & vField.mesh().Sf()
    )
    );
    Info << nl << "create dual flux field " << coarseField.name()+"_dual"<< endl;
    return phiFields_.back();
}


void Foam::dualMesh::update()
{
    for(label i=0;i<scalarFields_.first().size();i++) {
        mapFieldCoarseToFine(scalarFields_.first().at(i), scalarFields_.second().at(i));
    }
    for(label i=0;i<vectorFields_.first().size();i++) {
        mapFieldCoarseToFine(vectorFields_.first().at(i), vectorFields_.second().at(i));
    }
    forAll(phiFields_, fieldi)
    {
        volVectorField& vField = vectorFields_.second().at(fieldi);
        phiFields_.at(fieldi) = linearInterpolate(vField) & vField.mesh().Sf();
    }
}

// ************************************************************************* //
