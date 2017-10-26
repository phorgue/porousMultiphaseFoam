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

#include "krBrooksAndCorey3p.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace triRelativePermeabilityModels
{
defineTypeNameAndDebug(krBrooksAndCorey3p, 0);

addToRunTimeSelectionTable
(
    triRelativePermeabilityModel,
    krBrooksAndCorey3p,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triRelativePermeabilityModels::krBrooksAndCorey3p::krBrooksAndCorey3p
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& Sw,
    const volScalarField& So
)
    :
    triRelativePermeabilityModel(name, transportProperties, Sw, So),
    krBrooksAndCorey3pCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    n_
    (
        IOobject
        (
            "n",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        krBrooksAndCorey3pCoeffs_.lookupOrDefault<scalar>("n",0)
    ),
    Sar_
    (
        IOobject
        (
            "Sar",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        transportProperties.lookupOrDefault("Sar",dimensionedScalar("Sar",dimless,0))
    ),
    Sor_
    (
        IOobject
        (
            "Sor",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        transportProperties.lookupOrDefault("Sor",dimensionedScalar("Sor",dimless,0))
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
    kramax_
    (
        IOobject
        (
            "kramax",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw.mesh(),
        transportProperties.lookupOrDefault<scalar>("kramax",1.0)
    ),
    kromax_
    (
        IOobject
        (
            "kromax",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw.mesh(),
        transportProperties.lookupOrDefault<scalar>("kromax",1.0)
    ),
    krwmax_
    (
        IOobject
        (
            "krwmax",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sw.mesh(),
        transportProperties.lookupOrDefault<scalar>("krwmax",1.0)
    ),
    Sae_(( 1 - So_ - Sw_ ) / (1 - Sor_ - Swr_ - Sar_)),
    Soe_(( So_ - Sor_ ) / (1 - Sor_ - Swr_ - Sar_)),
    Swe_(( Sw_ - Swr_ ) / (1 - Sor_ - Swr_ - Sar_))
{
    if (gMin(n_) <= 0)
    {
        FatalErrorIn
            (
                "in krBrooksAndCorey3p.C"
            )
            << "3-phases relative permeability coefficient 'n' equal or less than 0" 
                << exit(FatalError);
    }
}

// ************************************************************************* //
