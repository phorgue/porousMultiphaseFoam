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
    Foam::relativePermeabilityModels::krIppisch

Description
    Ippisch Vogel Bastian relative permeability model (AWR 2006)

SourceFiles
    krIppisch.C

\*---------------------------------------------------------------------------*/

#ifndef krIppisch_H
#define krIppisch_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
    Class krIppisch Declaration
    \*---------------------------------------------------------------------------*/

class krIppisch
:
    public relativePermeabilityModel
{
    //- Private data
    volScalarField Smin_;
    volScalarField Smax_;

    //- Van Genuchten coefficients
    dictionary krIppischCoeffs_;
    volScalarField m_,n_,alpha_,tau_,he_;

    //- effective saturation of phase b
    volScalarField Se_;

    //- cut-off point saturation (computed with entry pressure)
    volScalarField Sc_;

public:

    //- Runtime type information
    TypeName("Ippisch");

    // Constructors

    //- Construct from components
    krIppisch
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sb
    );

    //- Destructor
    ~krIppisch()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct()
        {
            Se_==(Sb_-Smin_)/(Smax_-Smin_);
            krb_ = pow(Se_,tau_) * pow(
                (1-pow(1-pow(Sc_*Se_,1/m_),m_)) /
                (1-pow(1-pow(Sc_,1/m_),m_))
                ,2);
        }
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
