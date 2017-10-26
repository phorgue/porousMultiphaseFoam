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

#include "pcParkerAndLenhard.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace triCapillarityModels
{
defineTypeNameAndDebug(pcParkerAndLenhard, 0);

addToRunTimeSelectionTable
(
    triCapillarityModel,
    pcParkerAndLenhard,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triCapillarityModels::pcParkerAndLenhard::pcParkerAndLenhard
(
    const word& name,
    const dictionary& triCapillarityProperties,
    const volScalarField& Sw,
    const volScalarField& So
)
    :
    triCapillarityModel(name, triCapillarityProperties,Sw, So),
    pcParkerAndLenhardCoeffs_(triCapillarityProperties.subDict(typeName + "Coeffs")),
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
        Sw.mesh(),
        pcParkerAndLenhardCoeffs_.lookupOrDefault<scalar>("m",0)
    ),
    n_(1/(1-m_)),
    pc0_
    (
        IOobject
        (
            "pc0",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw.mesh(),
        pcParkerAndLenhardCoeffs_.lookupOrDefault("pc0",dimensionedScalar("pc0",dimensionSet(1,-1,-2,0,0),0.))
    ),
    beta_ao_
    (
        IOobject
        (
            "beta_ao",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw.mesh(),
        pcParkerAndLenhardCoeffs_.lookupOrDefault<scalar>("beta_ao",0)
    ),
    beta_ow_
    (
        IOobject
        (
            "beta_ow",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw.mesh(),
        pcParkerAndLenhardCoeffs_.lookupOrDefault<scalar>("beta_ow",0)
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
        pcParkerAndLenhardCoeffs_.lookupOrDefault("Swr",dimensionedScalar("Swr",dimless,0))
    ),
    Swe_((Sw_ - Swr_)/(1-Swr_)),
    St_((Sw_ + So_ - Swr_) / (1-Swr_))
{
    if (gMin(m_) == 0) FatalErrorIn("Foam::triCapillarityModels::pcParkerAndLenhard::pcParkerAndLenhard") << "m = 0 in pcParkerAndLenhard" << abort(FatalError);
    if (gMin(beta_ao_) <= 0)  FatalErrorIn("Foam::triCapillarityModels::pcParkerAndLenhard::pcParkerAndLenhard") << "beta_ao <= 0 in pcParkerAndLenhard" << abort(FatalError);    
    if (gMin(beta_ow_) <= 0)  FatalErrorIn("Foam::triCapillarityModels::pcParkerAndLenhard::pcParkerAndLenhard") << "beta_ow <= 0 in pcParkerAndLenhard" << abort(FatalError);
    
    correct();
    
}

// ************************************************************************* //
