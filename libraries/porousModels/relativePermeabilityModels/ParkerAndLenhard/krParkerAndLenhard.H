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
    Foam::triRelativePermeabilityModels::krParkerAndLenhard

Description
    Standard Parker and Lenhard (1990) three-phases relative permeability model.

SourceFiles
    krParkerAndLenhard.C

\*---------------------------------------------------------------------------*/

#ifndef krParkerAndLenhard_H
#define krParkerAndLenhard_H

#include "triRelativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace triRelativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
                       Class krParkerAndLenhard Declaration
\*---------------------------------------------------------------------------*/

class krParkerAndLenhard
:
    public triRelativePermeabilityModel
{
    // Private data
    dictionary krParkerAndLenhardCoeffs_;
    
    //--- coefficient
    volScalarField m_;
    
    //--- Residual water saturation
    volScalarField Swr_;

    //--- Effective and total saturations
    volScalarField Swe_;
    volScalarField St_;

public:

    //- Runtime type information
    TypeName("ParkerAndLenhard");

    // Constructors

    //- Construct from components
    krParkerAndLenhard
    (
        const word& name,
        const dictionary& triRelativePermeabilityProperties,
        const volScalarField& Sw,
        const volScalarField& So
    );

    //- Destructor
    ~krParkerAndLenhard()
        {}

    // Member Functions
    
    //- Correct the relative permeabilities
    void correct()
        {
            Swe_ == (Sw_ - Swr_)/(1-Swr_);
            St_ == (Sw_ + So_ - Swr_) / (1-Swr_);
            
            kra_ = pow(1-St_,0.5) * pow(1-pow(St_,1/m_),2*m_);
            kro_ = pow(St_ - Swe_,0.5) * pow( pow(1-pow(Swe_,1/m_),m_) - pow(1-pow(St_,1/m_),m_) ,2);
            krw_ = pow(Swe_,0.5) * pow(1-pow(1-pow(Swe_,1/m_),m_),2);
            
            dkradS_ = - pow((1-pow(St_,1/m_)),2*m_-1) * (-5*pow(St_,1/m_+1)+4*pow(St_,1/m_)+St_);
            dkradS_ *= 1/(2*pow((1-St_),0.5)*St_);
            dkradS_ *=  1/(1 - Swr_);

            dkrwdS_ = 0.5 * (1-pow((1-pow(Swe_,1/m_)),m_));
            dkrwdS_ *= ( 4 * pow(Swe_,1/m_-1/2) * pow( (1-pow(Swe_,1/m_)) , m_-1)) - (pow((1-pow(Swe_,1/m_)),m_)-1)/pow(Swe_,0.5);
            dkrwdS_ *= 1/(1 - Swr_);
        }
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace triRelativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
