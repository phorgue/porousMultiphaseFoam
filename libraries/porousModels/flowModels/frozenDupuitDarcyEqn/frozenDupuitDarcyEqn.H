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
    Foam::flowModels::frozenDupuitDarcyEqn

Description
    frozenDupuit-Darcy equation for 2D groundwater flow simulations

SourceFiles
    frozenDupuitDarcyEqn.C

\*---------------------------------------------------------------------------*/

#ifndef frozenDupuitDarcyEqn_H
#define frozenDupuitDarcyEqn_H

#include "porousMediumModel.H"
#include "fluidPhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam::flowModels
{

/*---------------------------------------------------------------------------*\
    Class frozenDupuitDarcyEqn Declaration
\*---------------------------------------------------------------------------*/

class frozenDupuitDarcyEqn
{

    //- Internal fields
    volScalarField hwater_;
    volScalarField seepageTerm_;

    //- Internal scalars
    const scalar zScale_;

    //- Internal surface fields
    surfaceScalarField phihwater_;

public:
    // Constructors

    //- Construct from components
    frozenDupuitDarcyEqn
    (
        const fvMesh& mesh,
        const IOdictionary& transportProperties,
        porousMediumModel& pmModel,
        fluidPhase& fluid
    );

    //- Destructor
    ~frozenDupuitDarcyEqn() = default;

    //- Access functions
    const volScalarField& hwater() const {return hwater_;}
    const surfaceScalarField& phihwater() const {return phihwater_;}
    const volScalarField& seepage() {return seepageTerm_;}
    const scalar& zScale() {return zScale_;}

    // Member Functions
    void info();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
