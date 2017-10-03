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

#include "alphaDispersion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
defineTypeNameAndDebug(alphaDispersion, 0);

addToRunTimeSelectionTable
(
    dispersionModel,
    alphaDispersion,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModels::alphaDispersion::alphaDispersion
(
    const word& name,
    const dictionary& transportProperties,
    const volVectorField& U
)
    :
    dispersionModel(name, transportProperties,U),
    alphaDispersionCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    tau_(alphaDispersionCoeffs_.lookup("tau")),
    alphaL_(alphaDispersionCoeffs_.lookup("alphaL")),
    alphaT_(alphaDispersionCoeffs_.lookup("alphaT"))
{
}

// ************************************************************************* //