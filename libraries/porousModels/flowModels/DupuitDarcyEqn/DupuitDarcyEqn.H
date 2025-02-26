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
    Foam::flowModels::DupuitDarcyEqn

Description
    Dupuit-Darcy equation for 2D groundwater flow simulations

SourceFiles
    DupuitDarcyEqn.C

\*---------------------------------------------------------------------------*/

#ifndef DupuitDarcyEqn_H
#define DupuitDarcyEqn_H

#include "fvMatrix.H"
#include "porousMediumModel.H"
#include "incompressiblePhase.H"
#include "multiDtManager.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam::flowModels
{

/*---------------------------------------------------------------------------*\
    Class DupuitDarcyEqn Declaration
\*---------------------------------------------------------------------------*/

class DupuitDarcyEqn
{
    //- gravity
    const dimensionedScalar g_;

    //- Internal fields
    volScalarField potential_;
    volScalarField z0_;
    volScalarField DEM_;
    volScalarField hwater_;
    volScalarField infiltration_;
    volScalarField seepageTerm_;

    //- References
    //--- mesh
    const fvMesh& mesh_;
    //--- porous medium models
    porousMediumModel& pmModel_;
    //--- water source term
    const volScalarField& sourceTerm_;
    //--- porosity/permeability field
    const volScalarField& eps_;
    const volScalarField& K_;
    //--- velocity field
    volVectorField& U_;

    //- Internal booleans
    const bool steady_, seepage_;

    //- Internal word
    const word DEMfileName_;

    //- Internal scalars
    const scalar zScale_;
    const dimensionedScalar rho_, mu_;
    const dimensionedScalar hwaterMin_;
    scalar cumulativeWaterAdded_, flowInOutFixedPoints_, flowOutSeepage_;

    //- Internal surface fields
    surfaceScalarField phi_, Kf_, Mf_, transmissivity_, phiG_, phiPc_, phihwater_;

    //- Internal variables for seepage/dry/fixed potential cell term
    labelList seepageIDList_, dryCellIDList_, fixedPotentialIDList_;
    scalarList seepageValueList_, fixedPotentialValueList_;
    volScalarField cellFlux_;

    //- Internal functions
    void initFixedPotential(const List<Tuple2<point,scalar>>&, bool);
    void correctInitialPotential();
    void updateProperties();
    void updateSeepage();
    void updateDryCells();

public:
    // Constructors

    //- Construct from components
    DupuitDarcyEqn
    (
        const fvMesh& mesh,
        const IOdictionary& transportProperties,
        porousMediumModel& pmModel,
        incompressiblePhase& fluidPhase,
        bool steady
    );

    //- Destructor
    ~DupuitDarcyEqn() = default;

    //- Access functions
    const volScalarField& potential() const {return potential_;}
    const volScalarField& hwater() const {return hwater_;}
    const surfaceScalarField& phi() const {return phi_;}
    const surfaceScalarField& phihwater() const {return phihwater_;}
    volScalarField& infiltration() {return infiltration_;}
    const volScalarField& seepage() {return seepageTerm_;}
    const labelList& dryCellIDList() const {return dryCellIDList_;}
    const scalar& zScale() {return zScale_;}
    const scalar& flowInOutFixedPoints() {return flowInOutFixedPoints_;}
    const scalar& flowOutSeepage() {return flowOutSeepage_;}

    // Member Functions
    scalar solve();
    void info();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
