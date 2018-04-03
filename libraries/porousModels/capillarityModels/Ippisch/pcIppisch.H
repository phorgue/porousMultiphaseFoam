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
    Foam::capillarityModels::Ippisch

Description
    Ippisch-Vogel-Bastian capillary pressure model (AWR 2006)

SourceFiles
    pcIppisch.C

\*---------------------------------------------------------------------------*/

#ifndef pcIppisch_H
#define pcIppisch_H

#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
                    Class pcIppisch Declaration
\*---------------------------------------------------------------------------*/

class pcIppisch
:
    public capillarityModel
{
    //- Ippisch coefficients
    dictionary pcIppischCoeffs_;
    volScalarField Smin_,Smax_;
    volScalarField m_,n_,alpha_,tau_,he_;

    // effective saturation of phase b
    volScalarField Se_;

    //- cut-off point saturation (computed with entry pressure)
    volScalarField Sc_;
    
public:

    //- Runtime type information
    TypeName("Ippisch");

    // Constructors

    //- Construct from components
    pcIppisch
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sb
    );

    //- Destructor
    ~pcIppisch()
        {}

    // Member Functions

    //- Correct the capillary pressure
    void correct()
        {
            FatalErrorIn("void correct()")
                << "Not implemented in pcIppish.H " << abort(FatalError);
        }

    //- Correct effective saturation and capillary capacity and return saturation (Richards' solver)
    volScalarField correctAndSb(volScalarField h)
        {
            //- Compute Saturation from head pressure
            h.dimensions().reset(dimless);
            Se_ = neg(h+he_) *(1/Sc_)*pow(1.0+pow(alpha_*mag(h),n_),-m_) + pos(h+he_+VSMALL) ;
            h.dimensions().reset(dimLength);

            //- Update capillary capacity
            Ch_.dimensions().reset(dimless);
            Ch_ == alpha_*m_*(Smax_-Smin_)/(1.0-m_)*pow(Se_,1.0/m_)*pow(1.0-pow(Se_,1.0/m_),m_) ;
            Ch_.dimensions().reset(dimless/dimLength);

            //- Update saturation
            volScalarField Stmp_ = Smin_+(Smax_-Smin_)*(Se_-SMALL);
    
            return Stmp_;
        }

    //- Capillary capacity (Richards' model)
    tmp<volScalarField> Ch() const
        {
            return Ch_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
