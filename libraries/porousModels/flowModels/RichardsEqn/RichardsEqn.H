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

Class
    Foam::flowModels::RichardsEqn

Description
    Richards equation for unsaturated groundwater water flows

SourceFiles
    RichardsEqn.C

\*---------------------------------------------------------------------------*/

#ifndef RichardsEqn_H
#define RichardsEqn_H

#include "fvMatrix.H"
#include "twophasePorousMediumModel.H"
#include "multiDtManager.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowModels
{

/*---------------------------------------------------------------------------*\
    Class RichardsEqn Declaration
\*---------------------------------------------------------------------------*/

class RichardsEqn
{
    //- gravity
    const uniformDimensionedVectorField g_;

    //- Internal fields
    volScalarField h_;
    volScalarField deltah_;
    volScalarField theta_;

    //- References
    //--- mesh
    const fvMesh& mesh_;
    //--- porous medium models
    twophasePorousMediumModel& pmModel_;
    relativePermeabilityModel& krModel_;
    capillarityModel& pcModel_;
    //--- water source term
    const volScalarField& sourceTerm_;
    //--- permeability field
    const volScalarField& K_;
    //--- velocity field
    volVectorField& U_;

    //- Internal boolean
    const bool massConservative_;
    const bool steady_;

    //- Internal scalars
    const dimensionedScalar rho_, mu_, Ss_;

    //- Internal surface fields
    surfaceScalarField phi_, Kf_, krf_, Lf_, Mf_, phiG_, phiPc_;

    //- Internal variables for seepage term
    const word patchDEM_;
    const label patchDEMID_;
    labelList topCellID_, seepageIDList_;
    scalarList distanceToDEM_, seepageValueList_;

    //- Internal functions
    fvScalarMatrix buildEqn();
    scalar initResidual(const fvScalarMatrix& hEqn);

public:
    // Constructors

    //- Construct from components
    RichardsEqn
    (
        const fvMesh& mesh,
        const IOdictionary& transportProperties,
        twophasePorousMediumModel& pmModel,
        incompressiblePhase& fluidPhase,
        const bool steady
    );

    //- Destructor
    virtual ~RichardsEqn()
        {}

    //- Access functions
    const volScalarField& h() {return h_;}
    volScalarField& theta() {return theta_;}
    surfaceScalarField& phi() {return phi_;}
    const labelList& seepageIDList() {return seepageIDList_;}

    // Member Functions
    void updateProperties(bool derivative);
    void updateSeepage();
    void noConvergence(multiDtManager& MDTM, Time& runTime, label algoID);
    const Tuple2<scalar, scalar> solvePicard(const scalar tolerance);
    const Tuple2<scalar, scalar> solveNewton(const scalar tolerance);
    void info();
    void checkSteadyConfig(const scalar tol1, const scalar tol2);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
