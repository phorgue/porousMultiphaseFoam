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
    Foam::capillarityModels::VanGenuchten

Description
    Standard Van Genuchten capillary pressure model

SourceFiles
    pcVanGenuchten.C

\*---------------------------------------------------------------------------*/

#ifndef pcVanGenuchten_H
#define pcVanGenuchten_H

#include "capillarityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  namespace capillarityModels
  {

    /*---------------------------------------------------------------------------*\
                        Class pcVanGenuchten Declaration
    \*---------------------------------------------------------------------------*/

    class pcVanGenuchten
      :
      public capillarityModel
    {
      // pcVan Genuchten coefficients
      dictionary pcVanGenuchtenCoeffs_;
      dimensionedScalar Sminpc_;    
      dimensionedScalar Smaxpc_;    
        scalar m_,n_;
      dimensionedScalar pc0_;

      // capillary pressure and derivative
      volScalarField pc_;
      volScalarField dpcdS_;
        
    public:

      //- Runtime type information
      TypeName("VanGenuchten");

      // Constructors

      //- Construct from components
      pcVanGenuchten
      (
       const word& name,
       const dictionary& capillarityProperties,
       const volScalarField& Sb
       );

      //- Destructor
      ~pcVanGenuchten()
      {}

      // Member Functions

      //- Return the capillary pressure according to pcVanGenuchten correlation
      tmp<volScalarField> pc() const
      {
	return pc_;
      }

      //- Return the capillary pressure derivative according to pcVanGenuchten correlation
      tmp<volScalarField> dpcdS() const
      {
	return dpcdS_;
      }

      //- Correct the capillarity
      void correct()
      {
	volScalarField Se_ = (Sb_- Sminpc_)/(Smaxpc_-Sminpc_);
        pc_ = pc0_ * pow(pow(Se_,-1/m_)-1,1/n_);
	dpcdS_ = - (1/(n_*m_)) * (pc0_/(Smaxpc_-Sminpc_)) * pow(pow(Se_,-1/m_)-1,(1/n_)-1) * pow(Se_,-(1+m_)/m_);
      }

      //- Read transportProperties dictionary
      bool read(const dictionary& capillarityProperties);
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  } // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
