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
      dimensionedScalar Sminpc_,Smaxpc_;    
      scalar m_,n_,alpha_;
      dimensionedScalar pc0_;

      // effective saturation of phase b 
      volScalarField Se_;

      // capillary pressure, derivative and capacity
        volScalarField pc_, dpcdS_, Ch_;
        
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

      //- Correct the capillary pressure
      void correct()
      {
	Se_ == (Sb_- Sminpc_)/(Smaxpc_-Sminpc_);
        pc_ == pc0_ * pow(pow(Se_,-1/m_)-1,1/n_);
	dpcdS_ = - (1/(n_*m_)) * (pc0_/(Smaxpc_-Sminpc_)) * pow(pow(Se_,-1/m_)-1,(1/n_)-1) * pow(Se_,-(1+m_)/m_);
      }

      //- Correct effective saturation and capillary capacity and return saturation (Richards' solver)
      volScalarField correctAndSb(volScalarField h) 
      {          
        //- Compute Saturation from! head pressure
        h.dimensions().reset(dimless);
        volScalarField Stmp_ = neg(h) *(Sminpc_+(Smaxpc_-Sminpc_)*pow(1.0+pow(alpha_*mag(h),n_),-m_)) + pos(h)*(Smaxpc_-SMALL);
        h.dimensions().reset(dimLength);  

        //- Update effective saturation and capillary capacity
	Se_ == (Stmp_- Sminpc_)/(Smaxpc_-Sminpc_);
        Ch_.dimensions().reset(dimless);
        Ch_ == alpha_*m_*(Smaxpc_-Sminpc_)/(1.0-m_)*pow(Se_,1.0/m_)*pow(1.0-pow(Se_,1.0/m_),m_) ;
        Ch_.dimensions().reset(dimless/dimLength);

        return Stmp_;
      }

      //- Capillary capacity (Richards' model)
      tmp<volScalarField> Ch() const
      {
        return Ch_;
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
