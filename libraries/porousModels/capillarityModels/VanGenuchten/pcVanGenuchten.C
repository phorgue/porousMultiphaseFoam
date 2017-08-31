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

#include "pcVanGenuchten.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{
defineTypeNameAndDebug(pcVanGenuchten, 0);

addToRunTimeSelectionTable
(
    capillarityModel,
    pcVanGenuchten,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcVanGenuchten::pcVanGenuchten
(
    const word& name,
    const dictionary& capillarityProperties,
    const volScalarField& Sb
)
    :
    capillarityModel(name, capillarityProperties,Sb),
    pcVanGenuchtenCoeffs_(capillarityProperties.subDict(typeName + "Coeffs")),
    Sminpc_
    (
        IOobject
        (
            Sb_.name()+"minpc",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        pcVanGenuchtenCoeffs_.lookupOrDefault(Sb_.name()+"minpc",capillarityProperties.lookupOrDefault(Sb_.name()+"min",dimensionedScalar(Sb_.name()+"min",dimless,0)))
    ),
    Smaxpc_
    (
        IOobject
        (
            Sb_.name()+"maxpc",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        pcVanGenuchtenCoeffs_.lookupOrDefault(Sb_.name()+"maxpc",capillarityProperties.lookupOrDefault(Sb_.name()+"max",dimensionedScalar(Sb_.name()+"max",dimless,0)))
    ),
    m_
    (
        IOobject
        (
            "m",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        pcVanGenuchtenCoeffs_.lookupOrDefault<scalar>("m",0)
    ),
    n_(1/(1-m_)),
    alpha_ // necessary for Richards solver
    (
        IOobject
        (
            "alpha",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        pcVanGenuchtenCoeffs_.lookupOrDefault<scalar>("alpha",0)
    ),
    pc0_
    (
        IOobject
        (
            "pc0",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        pcVanGenuchtenCoeffs_.lookupOrDefault("pc0",dimensionedScalar("pc0",dimensionSet(1,-1,-2,0,0),0.))
    ),
    Se_((Sb_- Sminpc_)/(Smaxpc_-Sminpc_))
{
    if (gMin(m_) == 0) FatalErrorIn("Foam::capillarityModels::pcVanGenuchten::pcVanGenuchten") << "m = 0 in pcVanGenuchten" << abort(FatalError);
    correct();
}

// ************************************************************************* //
