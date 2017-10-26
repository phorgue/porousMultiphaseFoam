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

#include "triCapillarityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(triCapillarityModel, 0);
defineRunTimeSelectionTable(triCapillarityModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triCapillarityModel::triCapillarityModel
(
    const word& name,
    const dictionary& triCapillarityProperties,
    const volScalarField& Sw,
    const volScalarField& So
)
    :
    name_(name),
    triCapillarityProperties_(triCapillarityProperties),
    Sw_(Sw),
    So_(So),
    pcao_
    (
        IOobject
        (
            "pcao",
            Sw.time().timeName(),
            Sw.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),       
        Sw.mesh(),
        dimensionSet(1,-1,-2,0,0,0,0)
    ),
    pcow_
    (
        IOobject
        (
            "pcow",
            Sw.time().timeName(),
            Sw.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),       
        Sw.mesh(),
        dimensionSet(1,-1,-2,0,0,0,0)
    ),
    pcaw_
    (
        IOobject
        (
            "pcaw",
            Sw.time().timeName(),
            Sw.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),       
        Sw.mesh(),
        dimensionSet(1,-1,-2,0,0,0,0)
    ),
    dpcaodS_
    (
        IOobject
        (
            "dpcaodS",
            Sw.time().timeName(),
            Sw.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),       
        Sw.mesh(),
        dimensionSet(1,-1,-2,0,0,0,0)
    ),    
    dpcowdS_
    (
        IOobject
        (
            "dpcowdS",
            Sw.time().timeName(),
            Sw.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),       
        Sw.mesh(),
        dimensionSet(1,-1,-2,0,0,0,0)
    ),    
    dpcawdS_
    (
        IOobject
        (
            "dpcawdS",
            Sw.time().timeName(),
            Sw.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),       
        Sw.mesh(),
        dimensionSet(1,-1,-2,0,0,0,0)
    )
{
}

// ************************************************************************* //
