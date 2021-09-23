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
    const word& name,
    const dictionary& capillarityProperties,
    const volScalarField& Sb
)
    :
    name_(name),
    capillarityProperties_(capillarityProperties),
    Sb_(Sb),
    pc_
    (
        IOobject
        (
            name+".pc",
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
            name+".dpcdS",
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
            name+".Ch",
            Sb.time().timeName(),
            Sb.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionSet(0,-1,0,0,0,0,0)
    ),
    Se_
    (
        IOobject
        (
            name+".Se",
            Sb.time().timeName(),
            Sb.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb,
        calculatedFvPatchScalarField::typeName
    )
{}

// ************************************************************************* //
