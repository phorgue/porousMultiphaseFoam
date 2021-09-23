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

    volScalarField Smin_;
    volScalarField Smax_;
    volScalarField pc0_;
    volScalarField pcMax_;

public:

    //- Runtime type information
    TypeName("linear");

    // Constructors

    //- Construct from components
    pcLinear
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sb
    );

    //- Destructor
    ~pcLinear()
        {}

    // Member Functions

    //- Correct the capillarity
    void correct()
        {
            Se_ = (Sb_- Smin_)/(Smax_-Smin_);
            pc_ = pc0_ + (1 - Se_) * (pcMax_ - pc0_);
            dpcdS_ =  pc0_ - pcMax_;
        }

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    volScalarField correctAndSb(volScalarField h)
        {
            FatalErrorIn("volScalarField correctAndSb(volScalarField h) const")
                << "Not implemented in pcLinear.H " << abort(FatalError);
            return Se_;
        };

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    scalar correctAndSb(const volScalarField& h, const label& celli)
        {
            FatalErrorIn("scalar correctOneCell(const volScalarField& h, const label& celli)")
                << "Not implemented in pcLinear.H " << abort(FatalError);
            return 0;
        };

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
