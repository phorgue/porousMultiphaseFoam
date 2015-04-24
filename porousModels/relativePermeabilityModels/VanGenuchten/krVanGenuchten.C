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

#include "krVanGenuchten.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace relativePermeabilityModels
  {
    defineTypeNameAndDebug(krVanGenuchten, 0);

    addToRunTimeSelectionTable
    (
     relativePermeabilityModel,
     krVanGenuchten,
     dictionary
     );
  }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krVanGenuchten::krVanGenuchten
(
 const word& name,
 const dictionary& relativePermeabilityProperties,
 const volScalarField& Sb
 )
  :
  relativePermeabilityModel(name, relativePermeabilityProperties,Sb),
  Smin_(relativePermeabilityProperties.lookup(Sb_.name()+"min")),
  Smax_(relativePermeabilityProperties.lookup(Sb_.name()+"max")),
  krVanGenuchtenCoeffs_(relativePermeabilityProperties.subDict(typeName + "Coeffs")),
  m_(krVanGenuchtenCoeffs_.lookupOrDefault<scalar>("m",0)),  
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
  kra_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
   Sb.mesh(),
   dimensionSet(0,0,0,0,0)
   ),
  krb_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb.mesh(),
   dimensionSet(0,0,0,0,0)
  ),
 dkradS_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
   Sb.mesh(),
   dimensionSet(0,0,0,0,0)
   ),
  dkrbdS_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb.mesh(),
   dimensionSet(0,0,0,0,0)
   ),
kramax_(krVanGenuchtenCoeffs_.lookupOrDefault<scalar>("kr"+Sb_.name()+"max",1.0)),
krbmax_(krVanGenuchtenCoeffs_.lookupOrDefault<scalar>("kr"+Sb_.name()+"max",1.0))
{
  if (m_ <= 0)
    {
      FatalErrorIn
        (
         "in krVanGenuchten.C"
         )
        << "Relative permeability coefficient m equal or less than 0" 
        << exit(FatalError);
    }
  correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::relativePermeabilityModels::krVanGenuchten::read
(
 const dictionary& relativePermeabilityProperties
 )
{
  relativePermeabilityProperties_ = relativePermeabilityProperties;

  Smin_=relativePermeabilityProperties.lookup(Sb_.name()+"min");
  Smax_=relativePermeabilityProperties.lookup(Sb_.name()+"max");

  krVanGenuchtenCoeffs_ = relativePermeabilityProperties.subDict(typeName + "Coeffs");
  m_ = krVanGenuchtenCoeffs_.lookupOrDefault<scalar>("m",0) ;

  if (m_ <= 0)
    {
      FatalErrorIn
	(
	 "in krVanGenuchten.C"
	 )
	<< "Van Genuchten coefficient m equal or less than 0" 
	<< exit(FatalError);
    }

  return true;
}


// ************************************************************************* //
