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

#include "porousMediumModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(porousMediumModel, 0);
defineRunTimeSelectionTable(porousMediumModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumModel::porousMediumModel
(
    const word Sname,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const autoPtr<incompressiblePhase>& phase,
    const word porousRegion
)
    :
    Sname_(Sname),
    transportProperties_(transportProperties),
    phase_(phase),
    eps_
    (
        IOobject
        (
            "eps",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        transportProperties.lookupOrDefault("eps",dimensionedScalar("",dimless,0.))
    ),
    K_
    (
        IOobject
        (
            "K",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    sourceTerm_
    (
        IOobject
        (
            "porousMediumSourceTerm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,0,-1,0,0),0.)
    )
{
    scalar Kfactor(transportProperties.getOrDefault<scalar>("Kfactor",1));
    if (Kfactor != 1)
    {
        K_ *= Kfactor;
        Info  << nl << "Reading permeability field factor : Kfactor = " << Kfactor << endl;
    }
    pcModel_ = capillarityModel::New(mesh, transportProperties, Sname, porousRegion);
    krModel_ = relativePermeabilityModel::New(mesh, transportProperties, Sname, porousRegion);
}

// ************************************************************************* //
