/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

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
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname
 )
    :
  capillarityModel(mesh, transportProperties, Sname),
  pcBrooksAndCoreyCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
  Smin_
  (
      IOobject
      (
          Sname+"min",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      mesh,
      transportProperties.getOrDefault(Sname+"min",dimensionedScalar(Sname+"min",dimless,0))
  ),
  Smax_
  (
      IOobject
      (
          Sname+"max",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      mesh,
      transportProperties.getOrDefault(Sname+"max",dimensionedScalar(Sname+"max",dimless,0))
  ),
  pc0_
  (
      IOobject
      (
          "pc0",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      mesh,
      pcBrooksAndCoreyCoeffs_.getOrDefault("pc0",dimensionedScalar("pc0",dimensionSet(1,-1,-2,0,0),0))
  ),
  hd_
  (
      IOobject
      (
          "hd",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("hd",dimLength,pcBrooksAndCoreyCoeffs_.getOrDefault<scalar>("hd",0))
  ),
  alpha_
  (
      IOobject
      (
          "alpha",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("alpha",dimless,pcBrooksAndCoreyCoeffs_.getOrDefault<scalar>("alpha",0))
  )
{
    if (gMin(alpha_) == 0) FatalErrorIn("Foam::capillarityModels::pcBrooksAndCorey::pcBrooksAndCorey") << "alpha = 0 in pcBrooksAndCorey" << abort(FatalError);

    Info << "Brooks and Corey parameters for capillary pressure model" << nl << "{" << endl;
    Info << "    pc0 ";
    if (pc0_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(pc0_).value() << endl;}
    Info << "    alpha ";
    if (alpha_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alpha_).value() << endl;}
    Info << "    hd ";
    if (hd_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(hd_).value() << endl;}
    Info <<  "    Smin ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    Smax ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info << "} \n" << endl;

}

// ************************************************************************* //
