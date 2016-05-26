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

#include "pcBrooksAndCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace capillarityModels
  {
    defineTypeNameAndDebug(pcBrooksAndCorey, 0);

    addToRunTimeSelectionTable
    (
     capillarityModel,
     pcBrooksAndCorey,
     dictionary
     );
  }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcBrooksAndCorey::pcBrooksAndCorey
(
 const word& name,
 const dictionary& capillarityProperties,
 const volScalarField& Sb
 )
    :
  capillarityModel(name, capillarityProperties,Sb),	
  pcBrooksAndCoreyCoeffs_(capillarityProperties.subDict(typeName + "Coeffs")),
  Sminpc_
  (
      IOobject
      (
          Sb_.name()+"minpc",
          "constant/porousModels",
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      pcBrooksAndCoreyCoeffs_.lookupOrDefault(Sb_.name()+"minpc",capillarityProperties.lookupOrDefault(Sb_.name()+"min",dimensionedScalar(Sb_.name()+"min",dimless,0)))
  ),
  Smaxpc_
  (
      IOobject
      (
          Sb_.name()+"maxpc",
          "constant/porousModels",
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
      ),
      Sb.mesh(),
      pcBrooksAndCoreyCoeffs_.lookupOrDefault(Sb_.name()+"maxpc",capillarityProperties.lookupOrDefault(Sb_.name()+"max",dimensionedScalar(Sb_.name()+"max",dimless,0)))
  ),
  pc0_
  (
      IOobject
      (
          "pc0",
          "constant/porousModels",
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      pcBrooksAndCoreyCoeffs_.lookupOrDefault("pc0",dimensionedScalar("pc0",dimless,0))
  ),
  alpha_
  (
      IOobject
      (
          "alpha",
          "constant/porousModels",
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      pcBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("alpha",0)
  ),
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
    if (gMin(alpha_) == 0) FatalErrorIn("Foam::capillarityModels::pcBrooksAndCorey::pcBrooksAndCorey") << "alpha = 0 in pcBrooksAndCorey" << abort(FatalError);
}

// ************************************************************************* //
