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

#include "krParkerAndLenhard.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace triRelativePermeabilityModels
{
defineTypeNameAndDebug(krParkerAndLenhard, 0);

addToRunTimeSelectionTable
(
    triRelativePermeabilityModel,
    krParkerAndLenhard,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triRelativePermeabilityModels::krParkerAndLenhard::krParkerAndLenhard
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& Sw,
    const volScalarField& So
)
    :
    triRelativePermeabilityModel(name, transportProperties, Sw, So),
    krParkerAndLenhardCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    m_
    (
        IOobject
        (
            "m",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        krParkerAndLenhardCoeffs_.lookupOrDefault<scalar>("m",0)
    ),
    Swr_
    (
        IOobject
        (
            "Swr",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        transportProperties.lookupOrDefault("Swr",dimensionedScalar("Swr",dimless,0))
    ),
    Swe_((Sw_ - Swr_)/(1-Swr_)),
    St_((Sw_ + So_ - Swr_) / (1-Swr_))
{
    if (gMin(m_) <= 0)
    {
        FatalErrorIn
            (
                "in krParkerAndLenhard.C"
            )
            << "triphase relative permeability coefficient m equal or less than 0" 
                << exit(FatalError);
    }
}

// ************************************************************************* //
