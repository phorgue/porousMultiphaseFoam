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

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(relativePermeabilityModel, 0);
defineRunTimeSelectionTable(relativePermeabilityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModel::relativePermeabilityModel
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname
)
    :
    Sname_(Sname),
    transportProperties_(transportProperties),
    kra_
    (
        IOobject
        (
            Sname+".kra",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0)
    ),
    krb_
    (
        IOobject
        (
            Sname+".krb",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0)
    ),
    dkradS_
    (
        IOobject
        (
            Sname+".dkradS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0)
    ),
    dkrbdS_
    (
        IOobject
        (
            Sname+".dkrbdS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0)
    ),
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
      dimensionedScalar(Sname+"min",dimless,0)
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
        dimensionedScalar(Sname+"min",dimless,1)
    ),
    Se_
    (
        IOobject
        (
            Sname+".Se",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0),
        calculatedFvPatchScalarField::typeName
    )
{}

// ************************************************************************* //
