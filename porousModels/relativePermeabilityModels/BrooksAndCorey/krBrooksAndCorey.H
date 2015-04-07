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
  Foam::relativePermeabilityModels::krBrooksAndCorey

  Description
  Standard krBrooks and Corey relativePermeability model.

  SourceFiles
  krBrooksAndCorey.C

  \*---------------------------------------------------------------------------*/

#ifndef krBrooksAndCorey_H
#define krBrooksAndCorey_H

#include "relativePermeabilityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  namespace relativePermeabilityModels
  {

    /*---------------------------------------------------------------------------*\
      Class krBrooksAndCorey Declaration
      \*---------------------------------------------------------------------------*/

    class krBrooksAndCorey
      :
      public relativePermeabilityModel
    {
      //- Private data
      dimensionedScalar Smin_;
      dimensionedScalar Smax_;

      dictionary krBrooksAndCoreyCoeffs_;
      scalar n_;

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
      TypeName("BrooksAndCorey");

      // Constructors

      //- Construct from components
      krBrooksAndCorey
      (
       const word& name,
       const dictionary& relativePermeabilityProperties,
       const volScalarField& Sb
       );

      //- Destructor
      ~krBrooksAndCorey()
      {}

      // Member Functions

      //- Return relative permeabilities according to the power law
      tmp<volScalarField> kra() const
      {
	return kra_;
      }   
      tmp<volScalarField> krb() const
      {
	return krb_;
      }
     
        //return the derivatives of relative permeabilities
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
	kra_ = kramax_* pow((scalar(1)-Se_),n_);
	krb_ = krbmax_ * pow(Se_,n_);


        dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
        dkrbdS_ = krbmax_*n_*pow(Se_,n_-1)/(Smax_- Smin_);	
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
