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
    Foam::triRelativePermeabilityModels::krBrooksAndCorey3p

Description
    "Modified" Brooks and Corey relative permeability model for 3-phase flow.

SourceFiles
    krBrooksAndCorey3p.C

\*---------------------------------------------------------------------------*/

#ifndef krBrooksAndCorey3p_H
#define krBrooksAndCorey3p_H

#include "triRelativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace triRelativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
                       Class krBrooksAndCorey3p Declaration
\*---------------------------------------------------------------------------*/

class krBrooksAndCorey3p
:
    public triRelativePermeabilityModel
{
    // Private data
    dictionary krBrooksAndCorey3pCoeffs_;

    //--- Power-law coefficient
    volScalarField n_;
    
    //--- Residual saturations
    volScalarField Sar_;
    volScalarField Sor_;
    volScalarField Swr_;
    
    //--- Kr max coefficients 
    volScalarField kramax_;
    volScalarField kromax_;    
    volScalarField krwmax_;
    
    //--- Effective saturations
    volScalarField Sae_;
    volScalarField Soe_;    
    volScalarField Swe_;

    
public:

    //- Runtime type information
    TypeName("BrooksAndCorey3p");

    // Constructors

    //- Construct from components
    krBrooksAndCorey3p
    (
        const word& name,
        const dictionary& triRelativePermeabilityProperties,
        const volScalarField& Sw,
        const volScalarField& So
    );

    //- Destructor
    ~krBrooksAndCorey3p()
        {}

    // Member Functions
    
    //- Correct the relative permeabilities
    void correct()
        {
            volScalarField Sae_ = ( 1 - So_ - Sw_ ) / (1 - Sor_ - Swr_ - Sar_);
            kra_ = kramax_ * pow(Sae_,n_);

            volScalarField Soe_ = ( So_ - Sor_ ) / (1 - Sor_ - Swr_ - Sar_);
            kro_ = kromax_ * pow(Soe_,n_);

            volScalarField Swe_ = ( Sw_ - Swr_ ) / (1 - Sor_ - Swr_ - Sar_);
            krw_ = krwmax_ * pow(Swe_,n_);        

            dkradS_ = kramax_*n_*pow(Soe_,n_-1)/(1 - Sor_ - Swr_ - Sar_);
            dkrwdS_ = krwmax_*n_*pow(Swe_,n_-1)/(1 - Sor_ - Swr_ - Sar_);  

        }
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace triRelativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
