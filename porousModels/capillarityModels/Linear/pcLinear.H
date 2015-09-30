/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
   \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

      dimensionedScalar Sminpc_;    
      dimensionedScalar Smaxpc_;    
      dimensionedScalar pc0_;
      dimensionedScalar pcMax_;

      //- effective saturation of phase b
      volScalarField Se_;

      // capillary pressure and derivative
      volScalarField pc_;
      volScalarField dpcdS_;
        
    public:

      //- Runtime type information
      TypeName("linear");

      // Constructors

      //- Construct from components
      pcLinear
      (
       const word& name,
       const dictionary& capillarityProperties,
       const volScalarField& Sb
       );

      //- Destructor
      ~pcLinear()
      {}

      // Member Functions

      //- Return the capillary pressure according to pcLinear correlation
      tmp<volScalarField> pc() const
      {
	return pc_;
      }

      //- Return the capillary pressure derivative according to pcLinear correlation
      tmp<volScalarField> dpcdS() const
      {
	return dpcdS_;
      }

      //- Correct the capillarity
      void correct()
      {
	Se_ == (Sb_- Sminpc_)/(Smaxpc_-Sminpc_);
        pc_ = pc0_ + (1 - Se_) * (pcMax_ - pc0_) ; 
	dpcdS_ =  pc0_ - pcMax_; 
      }

      //- Compute Moisture from head pressure (Richards' model)
      volScalarField correctAndSb(volScalarField h) 
      {
        FatalErrorIn("volScalarField correctAndSb(volScalarField h) const")
                   << "Not implemented in pcLinear.H " << abort(FatalError);
        return Se_;
      };

      //- Capillary capacity (Richards' model)
      tmp<volScalarField> Ch() const
      {
        FatalErrorIn("tmp<volScalarField> Ch() const")
                  << "Not implemented in pcLinear.H " << abort(FatalError);
        return Se_;
      };

      //- Read transportProperties dictionary
      bool read(const dictionary& capillarityProperties);

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  } // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
