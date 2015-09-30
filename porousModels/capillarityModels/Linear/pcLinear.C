/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
 const dictionary& capillarityProperties,
 const volScalarField& Sb
 )
    :
  capillarityModel(name, capillarityProperties,Sb),
  pcLinearCoeffs_(capillarityProperties.subDict(typeName + "Coeffs")),
  Sminpc_(pcLinearCoeffs_.lookupOrDefault(Sb_.name()+"minpc",dimensionedScalar(Sb_.name()+"min",capillarityProperties.lookup(Sb_.name()+"min")+0))),
  Smaxpc_(pcLinearCoeffs_.lookupOrDefault(Sb_.name()+"maxpc",dimensionedScalar(Sb_.name()+"max",capillarityProperties.lookup(Sb_.name()+"max")+0))),
  pc0_(pcLinearCoeffs_.lookup("pc0")),
  pcMax_(pcLinearCoeffs_.lookup("pcMax")),
  Se_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb_
   ),
  pc_
  (
   IOobject
   (
    name,
    Sb.time().timeName(),
    Sb.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb.mesh(),
   dimensionSet(1,-1,-2,0,0,0,0)
  ),
  dpcdS_
  (
   IOobject
   (
    name,
    Sb.time().timeName(),
    Sb.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb.mesh(),
   dimensionSet(1,-1,-2,0,0,0,0)
   )
{
  correct();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::capillarityModels::pcLinear::read
(
 const dictionary& capillarityProperties
 )
{
  capillarityProperties_ = capillarityProperties;

  pcLinearCoeffs_ = capillarityProperties.subDict(typeName + "Coeffs");
  pcLinearCoeffs_.lookup(Sb_.name()+"minpc") >> Sminpc_;
  pcLinearCoeffs_.lookup(Sb_.name()+"maxpc") >> Smaxpc_;
  pcLinearCoeffs_.lookup("pc0") >> pc0_;
  pcLinearCoeffs_.lookup("pcMax") >> pcMax_;
  return true;
}

// ************************************************************************* //
