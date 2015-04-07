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
  Foam::relativePermeabilityModels::krVanGenuchten

  Description
  Standard Brooks and Corey relativePermeability model.

  SourceFiles
  krVanGenuchten.C

  \*---------------------------------------------------------------------------*/

#ifndef krVanGenuchten_H
#define krVanGenuchten_H

#include "relativePermeabilityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  namespace relativePermeabilityModels
  {

    /*---------------------------------------------------------------------------*\
      Class krVanGenuchten Declaration
      \*---------------------------------------------------------------------------*/

    class krVanGenuchten
      :
      public relativePermeabilityModel
    {
      //- Private data
      dimensionedScalar Smin_;
      dimensionedScalar Smax_;

      //- Van Genuchten coefficients
      dictionary krVanGenuchtenCoeffs_;
      scalar m_;

      //- effective saturation of phase b
      volScalarField Se_;

      //- relative permeability of each phase 
      volScalarField kra_;
      volScalarField krb_;

      volScalarField dkradS_;
      volScalarField dkrbdS_;

      //- end points
      scalar kramax_;
      scalar krbmax_;


    public:

      //- Runtime type information
      TypeName("VanGenuchten");

      // Constructors

      //- Construct from components
      krVanGenuchten
      (
       const word& name,
       const dictionary& relativePermeabilityProperties,
       const volScalarField& Sb
       );

      //- Destructor
      ~krVanGenuchten()
      {}

      // Member Functions

      //- Return relative permeabilities according to the Van Genuchten law
      tmp<volScalarField> kra() const
      {
	return kra_;
      }   

      tmp<volScalarField> krb() const
      { 
	return krb_;
      }
       
      //- Return relative permeabilities according to the Van Genuchten law
      tmp<volScalarField> dkradS() const
      {
	return dkradS_;
      }   

      tmp<volScalarField> dkrbdS() const
      { 
	return dkrbdS_;
      }
      
 
      //- Correct the relative permeabilities
      void correct()
      {
	Se_==(Sb_-Smin_)/(Smax_-Smin_);
	kra_ = kramax_ * pow(1-Se_,0.5) * pow(1-pow(Se_,1/m_),2*m_);
	krb_ = krbmax_ * pow(Se_,0.5) * pow(1-pow(1-pow(Se_,1/m_),m_),2);
      
	dkradS_ = - pow((1-pow(Se_,1/m_)),2*m_-1) * (-5*pow(Se_,1/m_+1)+4*pow(Se_,1/m_)+Se_);
	dkradS_ *= 1/(2*pow((1-Se_),0.5)*Se_);
	dkradS_ *=  1/(Smax_ - Smin_);
	dkradS_ *= kramax_;

	 dkrbdS_ = 0.5 * (1-pow((1-pow(Se_,1/m_)),m_));
	 dkrbdS_ *= ( 4 * pow(Se_,1/m_-1/2) * pow( (1-pow(Se_,1/m_)) , m_-1)) - ( pow((1-pow(Se_,1/m_)),m_) -1) / pow(Se_,0.5);
	 dkrbdS_ *= 1/(Smax_ - Smin_);
	 dkrbdS_ *= krbmax_;

      }

      //- Read transportProperties dictionary
      bool read(const dictionary& relativePermeabilityProperties);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  } // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
