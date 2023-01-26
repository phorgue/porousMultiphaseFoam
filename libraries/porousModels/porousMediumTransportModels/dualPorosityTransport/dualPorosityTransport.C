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
    porousMediumTransportModel,
    dualPorosityTransport,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumTransportModels::dualPorosityTransport::dualPorosityTransport
(
    const porousMediumModel& pmModel
)
    :
    porousMediumTransportModel(pmModel),
    composition_(
        transportProperties_,
        wordList(transportProperties_.lookupOrDefault("species", wordList(1, "C"))),
        pmModel.mesh(),
        word::null,
        pmModel.eps(),
        &sourceEventList_,
        "C",
        dimless,
        "Matrix"
    ),
    UMatrix_(transportProperties_.db().lookupObject<volVectorField>("UMatrix")),
    phiMatrix_(transportProperties_.db().lookupObject<surfaceScalarField>("phiMatrix")),
    thetaMatrix_(transportProperties_.db().lookupObject<volScalarField>("thetaMatrix"))
{
}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

void Foam::porousMediumTransportModels::dualPorosityTransport::correct
(
)
{
    //- Correct matrix dispersion
    composition_.correct<volScalarField>(UMatrix_, thetaMatrix_);

    forAll(composition_.Y(), speciesi)
    {
        const auto& speciesName = composition_.species()[speciesi];
    
        auto& C = composition_.Y(speciesi);
        const auto& R = composition_.R(speciesi);
        const auto& Deff = composition_.Deff(speciesi);
        const auto& lambda = composition_.lambda(speciesi);
        const auto& sourceTerm = composition_.sourceTerm(speciesi);

        fvScalarMatrix CMatrixEqn
            (
                R * fvm::ddt(thetaMatrix_, C)
                + fvm::div(phiMatrix_, C, "div(phi,C)")
                -  fvm::laplacian(thetaMatrix_* Deff, C, "laplacian(Deff,C)")
                ==
                - sourceTerm
                - thetaMatrix_ * R * fvm::Sp(lambda,C)
            );

        CMatrixEqn.solve(pmModel_.mesh().solver("C"));

        Info<< "Concentration: Min(" << speciesName << ") = " << gMin(C.internalField())
            << " Max(" << speciesName << ") = " << gMax(C.internalField())
            << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*thetaMatrix_).value()
            << endl;
    }
}

// ************************************************************************* //
