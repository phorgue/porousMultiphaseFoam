/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "incompressiblePhase.H"
#include "fixedValueFvPatchFields.H"
#include "linear.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressiblePhase::incompressiblePhase
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
:  
    name_(phaseName),
    dict_(transportProperties.subDict("phase." + phaseName)),
    mu_(dict_.lookup("mu")),
    rho_(dict_.lookup("rho")),
    U_
    (
        IOobject
        (
            "U" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{   
    const word phiName = "phi" + phaseName;

    IOobject phiHeader
    (
        phiName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (phiHeader.headerOk())
    {
        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        for (label i=0; i<U_.boundaryField().size(); i++)
        {
            if (isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i]))
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Foam::linearInterpolate(U_) & mesh.Sf(),
                phiTypes
            )
        );
    }
}


Foam::autoPtr<Foam::incompressiblePhase> Foam::incompressiblePhase::New
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
{
    return autoPtr<incompressiblePhase>
    (
        new incompressiblePhase(mesh, transportProperties, phaseName)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressiblePhase::~incompressiblePhase()
{}


// ************************************************************************* //
