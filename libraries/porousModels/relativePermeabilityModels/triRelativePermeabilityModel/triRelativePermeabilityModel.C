/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "triRelativePermeabilityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(triRelativePermeabilityModel, 0);
defineRunTimeSelectionTable(triRelativePermeabilityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triRelativePermeabilityModel::triRelativePermeabilityModel
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& Sw,
    const volScalarField& So
)
    :
    name_(name),
    transportProperties_(transportProperties),
    Sw_(Sw),
    So_(So),
    kra_
    (
        IOobject
        (
            "kra",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Sw_.mesh(),
        dimensionSet(0,0,0,0,0)
    ),
    kro_
    (
        IOobject
        (
            "kro",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Sw_.mesh(),
        dimensionSet(0,0,0,0,0)
    ),
    krw_
    (
        IOobject
        (
            "krw",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Sw_.mesh(),
        dimensionSet(0,0,0,0,0)
    ),
    dkradS_
    (
        IOobject
        (
            "dkradS",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        dimensionSet(0,0,0,0,0)
    ),
    dkrwdS_
    (
        IOobject
        (
            "dkrwdS",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        dimensionSet(0,0,0,0,0)
    )  
{}

// ************************************************************************* //
