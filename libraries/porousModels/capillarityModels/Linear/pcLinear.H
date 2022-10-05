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
    Foam::capillarityModels::linear

Description
    Standard linear model for capillary pressure

SourceFiles
    pcLinear.C

\*---------------------------------------------------------------------------*/

#ifndef pcLinear_H
#define pcLinear_H

#include "capillarityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
    Class pcLinear Declaration
    \*---------------------------------------------------------------------------*/

class pcLinear
:
    public capillarityModel
{
    // pcVan Genuchten coefficients
    dictionary pcLinearCoeffs_;
    volScalarField pc0_;
    volScalarField pcMax_;

public:

    //- Runtime type information
    TypeName("linear");

    // Constructors

    //- Construct from components
    pcLinear
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname
    );

    //- Destructor
    ~pcLinear()
        {}

    // Member Functions

    //- Correct the capillarity
    void correct(const volScalarField& Sb)
        {
            Se_ = (Sb - Smin_)/(Smax_-Smin_);
            pc_ = pc0_ + (1 - Se_) * (pcMax_ - pc0_);
            dpcdS_ =  pc0_ - pcMax_;
        }

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    tmp<volScalarField> correctAndSb(volScalarField& h)
        {
            FatalErrorIn("volScalarField correctAndSb(volScalarField h) const")
                << "Not implemented in pcLinear.H " << abort(FatalError);
            return NULL;
        };

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    scalar correctAndSb(volScalarField& h, const label& celli)
        {
            FatalErrorIn("scalar correctAndSb(const volScalarField& h, const label& celli)")
                << "Not implemented in pcLinear.H " << abort(FatalError);
            return 0;
        };

    //- Compute saturation from given head pressure
    tmp<volScalarField> S(volScalarField& h)
        {
            FatalErrorIn("tmp<volScalarField> S(const volScalarField& h)")
                << "Not implemented in pcLinear.H " << abort(FatalError);
            return 0;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
