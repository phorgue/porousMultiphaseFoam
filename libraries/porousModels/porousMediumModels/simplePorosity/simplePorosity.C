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

#include "simplePorosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumModels
{
defineTypeNameAndDebug(simplePorosity, 0);

addToRunTimeSelectionTable
(
    porousMediumModel,
    simplePorosity,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumModels::simplePorosity::simplePorosity
(
    const word Sname,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const autoPtr<incompressiblePhase>& phase
)
    :
    porousMediumModel(Sname, mesh, transportProperties, phase)
{}

void Foam::porousMediumModels::simplePorosity::correct()
{
    //- nothing for simple porosity model
}

void Foam::porousMediumModels::simplePorosity::correct(const volScalarField& h, const bool steady, const bool conservative)
{
    //- nothing for simple porosity model
}
// ************************************************************************* //
