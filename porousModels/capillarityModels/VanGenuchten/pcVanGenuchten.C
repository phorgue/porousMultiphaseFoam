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
  Sminpc_(pcVanGenuchtenCoeffs_.lookupOrDefault(Sb_.name()+"minpc",dimensionedScalar(Sb_.name()+"min",capillarityProperties.lookup(Sb_.name()+"min")+0))),
  Smaxpc_(pcVanGenuchtenCoeffs_.lookupOrDefault(Sb_.name()+"maxpc",dimensionedScalar(Sb_.name()+"max",capillarityProperties.lookup(Sb_.name()+"max")+0))),
  m_(pcVanGenuchtenCoeffs_.lookupOrDefault<scalar>("m",0)),
  n_(1/(1-m_)),
  alpha_(pcVanGenuchtenCoeffs_.lookupOrDefault<scalar>("alpha",0)), // necessary for Richards solver
  pc0_(pcVanGenuchtenCoeffs_.lookupOrDefault("pc0",dimensionedScalar("pc0",dimensionSet(1,-1,-2,0,0),0))),
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
  ),
  Ch_
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
   dimensionSet(0,-1,0,0,0,0,0)
  )
{
        Info << " Saturation min = " << Sminpc_.value()
         << nl << " Saturation max = " << Smaxpc_.value() 
         << nl << " m = " << m_
         << nl << " pc0 = " << pc0_.value()
         << nl << " alpha = " << alpha_ << nl <<  endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::capillarityModels::pcVanGenuchten::read
(
 const dictionary& capillarityProperties
 )
{
  capillarityProperties_ = capillarityProperties;

  pcVanGenuchtenCoeffs_ = capillarityProperties.subDict(typeName + "Coeffs");
  pcVanGenuchtenCoeffs_.lookup(Sb_.name()+"minpc") >> Sminpc_;
  pcVanGenuchtenCoeffs_.lookup(Sb_.name()+"maxpc") >> Smaxpc_;
  pcVanGenuchtenCoeffs_.lookup("pc0") >> pc0_;
  m_ = pcVanGenuchtenCoeffs_.lookupOrDefault<scalar>("m",0);
  n_ = pcVanGenuchtenCoeffs_.lookupOrDefault<scalar>("n",0);
  return true;
}

// ************************************************************************* //
