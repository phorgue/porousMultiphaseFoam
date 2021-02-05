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
    Foam::capillarityModels::pcBrooksAndCorey

Description
    Standard Brooks and Corey capillary pressure model.

SourceFiles
    pcBrooksAndCorey.C

\*---------------------------------------------------------------------------*/

#ifndef pcBrooksAndCorey_H
#define pcBrooksAndCorey_H

#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
    Class pcBrooksAndCorey Declaration
    \*---------------------------------------------------------------------------*/

class pcBrooksAndCorey
:
    public capillarityModel
{
    // pcBrooks and Corey coefficients
    dictionary pcBrooksAndCoreyCoeffs_;
       
    volScalarField Smin_;
    volScalarField Smax_;
    volScalarField pc0_;
    volScalarField hd_;
    volScalarField alpha_;
        
public:

    //- Runtime type information
    TypeName("BrooksAndCorey");

    // Constructors

    //- Construct from components
    pcBrooksAndCorey
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sb
    );

    //- Destructor
    ~pcBrooksAndCorey()
        {}

    // Member Functions

    //- Correct the capillary pressure
    void correct()
        {
            Se_ = (Sb_-Smin_)/(Smax_ - Smin_);
            pc_ = pc0_ * pow(Se_,-alpha_);
            dpcdS_= -alpha_*pc0_*(pow(Se_+VSMALL,-alpha_-1))/(Smax_-Smin_);
        }

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    volScalarField correctAndSb(volScalarField h) 
        {
            if (gMax(hd_) == 0)
            {
                FatalErrorIn("pcBrooksAndCorey.H")
                    << "hd equal to zero in Brooks and Corey capillary pressure model"
                        << abort(FatalError);
            }

            //- Compute saturation from! head pressure
            h.dimensions().reset(dimless);
            Se_ = neg(h-hd_) * pow(hd_/(h+VSMALL),alpha_) + pos(h-hd_+VSMALL);
            h.dimensions().reset(dimLength);

            //- Update  capillary capacity
            Ch_.dimensions().reset(dimless);
            Ch_ = neg(h-hd_) * (alpha_ * (Smax_-Smin_) / mag(hd_)) * pow(hd_/(h+VSMALL),alpha_+1) + pos(h-hd_+VSMALL);
            Ch_.dimensions().reset(dimless/dimLength);

            //- Update saturation
            volScalarField Stmp_(Smin_+(Smax_-Smin_)*Se_);
            return Se_;
        };

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    scalar correctAndSb(const volScalarField& h, const label& celli)
        {
            //- Compute saturation from head pressure
            scalar Se = neg(h[celli]-hd_[celli]) * pow(hd_[celli]/(h[celli]+VSMALL),alpha_[celli]) + pos(h[celli]-hd_[celli]+VSMALL);

            //- Update capillary capacity
            Ch_[celli] = neg(h[celli]-hd_[celli]) * (alpha_[celli] * (Smax_[celli]-Smin_[celli]) / mag(hd_[celli])) * pow(hd_[celli]/(h[celli]+VSMALL),alpha_[celli]+1) + pos(h[celli]-hd_[celli]+VSMALL);

            //- Return moisture value
            return Smin_[celli]+(Smax_[celli]-Smin_[celli])*Se;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
