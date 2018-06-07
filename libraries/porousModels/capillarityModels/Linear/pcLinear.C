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

#include "pcLinear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{
defineTypeNameAndDebug(pcLinear, 0);

addToRunTimeSelectionTable
(
    capillarityModel,
    pcLinear,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcLinear::pcLinear
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& Sb
)
    :
    capillarityModel(name, transportProperties,Sb),
    pcLinearCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    Smin_
    (
        IOobject
        (
            Sb_.name()+"min",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        transportProperties.lookupOrDefault(Sb_.name()+"min",dimensionedScalar(Sb_.name()+"min",dimless,0))
    ),
    Smax_
    (
        IOobject
        (
            Sb_.name()+"max",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
       transportProperties.lookupOrDefault(Sb_.name()+"max",dimensionedScalar(Sb_.name()+"max",dimless,0))
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
        pcLinearCoeffs_.lookupOrDefault("pc0",dimensionedScalar("pc0",dimensionSet(1,-1,-2,0,0),0))
    ),
    pcMax_
    (
        IOobject
        (
            "pcMax",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        pcLinearCoeffs_.lookupOrDefault("pcMax",dimensionedScalar("pcMax",dimensionSet(1,-1,-2,0,0),0))
    ),
    Se_((Sb_- Smin_)/(Smax_-Smin_))
{

    Info << "Linear parameters for capillary pressure model" << nl << "{" << endl;
    Info << "    pc0 ";
    if (pc0_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(pc0_).value() << endl;}
    Info << "    pcMax ";
    if (pcMax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(pcMax_).value() << endl;}
    Info <<  "    Smin ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    Smax ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info << "} \n" << endl;
    
}

// ************************************************************************* //
