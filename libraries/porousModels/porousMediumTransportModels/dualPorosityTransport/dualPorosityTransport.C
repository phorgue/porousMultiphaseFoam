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

#include "fvCFD.H"
#include "dualPorosityTransport.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumTransportModels
{
defineTypeNameAndDebug(dualPorosityTransport, 0);

addToRunTimeSelectionTable
(
    porousMediumModel,
    dualPorosityTransport,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumTransportModels::dualPorosityTransport::dualPorosityTransport
(
    const word Sname,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const porousMediumModel& pmModel
)
    :
    porousMediumModelTransport(Sname, mesh, transportProperties, pmModel),
    multiscalarMixture(
        transportProperties,
        wordList(transportProperties.lookupOrDefault("species", wordList(1, "C"))),
        mesh,
        word::null,
        eps,
        &sourceEventList,
        "C",
        dimless,
        "Matrix"
    );
{
}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

void Foam::porousMediumModels::dualPorosityTransport::correct()
{
  
    // forAll(composition.Y(), speciesi)
    // {
    //     const auto& speciesName = composition.species()[speciesi];
    
    //     auto& C = composition.Y(speciesi);
    //     const auto& R = composition.R(speciesi);
    //     const auto& Deff = composition.Deff(speciesi);
    //     const auto& lambda = composition.lambda(speciesi);
    //     const auto& sourceTerm = composition.sourceTerm(speciesi);

    //     fvScalarMatrix CEqn
    //         (
    //             eps * R * Saturation * fvm::ddt(C)
    //             + fvm::div(phi, C, "div(phi,C)")
    //             -  fvm::laplacian(eps * Saturation * Deff, C, "laplacian(Deff,C)")
    //             ==
    //             - sourceTerm
    //             - eps * R * Saturation * fvm::Sp(lambda,C)
    //         );

    //     CEqn.solve(mesh.solver("C"));

    //     dtManager[speciesi].updateDerivatives();

    //     Info<< "Concentration: Min(" << speciesName << ") = " << gMin(C.internalField()) 
    //         << " Max(" << speciesName << ") = " << gMax(C.internalField())
    //         << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*Saturation*eps).value()
    //         << " dCmax = " << dtManager[speciesi].dVmax()*runTime.deltaTValue()
    //         << endl;

    // }
}

// ************************************************************************* //
