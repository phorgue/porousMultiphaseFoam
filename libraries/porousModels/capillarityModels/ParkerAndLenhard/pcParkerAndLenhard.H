/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::triCapillarityModels::ParkerAndLenhard

Description
    Standard Parker and Lenhard capillary pressure model

SourceFiles
    pcParkerAndLenhard.C

\*---------------------------------------------------------------------------*/

#ifndef pcParkerAndLenhard_H
#define pcParkerAndLenhard_H

#include "triCapillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace triCapillarityModels
{

/*---------------------------------------------------------------------------*\
    Class pcParkerAndLenhard Declaration
    \*---------------------------------------------------------------------------*/

class pcParkerAndLenhard
:
    public triCapillarityModel
{
    //- Private data
    dictionary pcParkerAndLenhardCoeffs_;
        
    //--- Parker and Lenhard coefficients
    volScalarField m_,n_;
    volScalarField pc0_,beta_ao_,beta_ow_;
      
    //--- Residual water saturation
    volScalarField Swr_;

    //--- effective water saturation and  total liquid saturation
    volScalarField Swe_,St_;
        
public:

    //- Runtime type information
    TypeName("ParkerAndLenhard");

    // Constructors

    //- Construct from components
    pcParkerAndLenhard
    (
        const word& name,
        const dictionary& triCapillarityProperties,
        const volScalarField& Sw,
        const volScalarField& So
    );

    //- Destructor
    ~pcParkerAndLenhard()
        {}

    // Member Functions
        
    //- Correct the capillary pressure
    void correct()
        {
            Swe_ == (Sw_ - Swr_)/(1-Swr_);
            St_ == (Sw_ + So_ - Swr_) / (1-Swr_);
          
            pcao_ = (pc0_/beta_ao_) * pow(pow(St_,-1/m_)-1,1/n_);
            pcow_ = (pc0_/beta_ow_) * pow(pow(Swe_,-1/m_)-1,1/n_);
            pcow_ = (pc0_) * pow(pow(Swe_,-1/m_)-1,1/n_);

            dpcaodS_ = - (pc0_/(n_*m_*(1-Swr_)*beta_ao_)) * pow(pow(St_,-1/m_)-1,(1/n_)-1) * pow(St_,-(1+m_)/m_);
            dpcowdS_ = - (pc0_/(n_*m_*(1-Swr_)*beta_ow_)) * pow(pow(Swe_,-1/m_)-1,(1/n_)-1) * pow(Swe_,-(1+m_)/m_);
            dpcawdS_ = - (pc0_/(n_*m_*(1-Swr_))) * pow(pow(Swe_,-1/m_)-1,(1/n_)-1) * pow(Swe_,-(1+m_)/m_);
        }

    //- Correct effective saturation and capillary capacity and return saturation (Richards' solver)
    // volScalarField correctAndSb(volScalarField h) 
    // { 
    // }

    //- Capillary capacity (Richards' model)
    // tmp<volScalarField> Ch() const
    // {
    //   return Ch_;
    // }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
