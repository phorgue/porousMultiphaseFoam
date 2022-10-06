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

#include "capillarityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(capillarityModel, 0);
defineRunTimeSelectionTable(capillarityModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModel::capillarityModel
(
    const fvMesh& mesh,
    const dictionary& capillarityProperties,
    const word& Sname,
    const word mediumName
)
    :
    Sname_(Sname),
    capillarityProperties_(capillarityProperties),
    pc_
    (
        IOobject
        (
            Sname+mediumName+".pc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1,-1,-2,0,0,0,0)
    ),
    dpcdS_
    (
        IOobject
        (
            Sname+mediumName+".dpcdS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1,-1,-2,0,0,0,0)
    ),
    Ch_
    (
        IOobject
        (
            Sname+mediumName+".Ch",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,-1,0,0,0,0,0)
    ),
    Smin_
    (
      IOobject
      (
          Sname+"min"+mediumName,
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar(dimless,0)
    ),
    Smax_
    (
        IOobject
        (
            Sname+"max"+mediumName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,1)
    ),
    Se_
    (
        IOobject
        (
            Sname+".Se"+mediumName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        calculatedFvPatchScalarField::typeName
    )
{}

// ************************************************************************* //
