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
    volScalarField alpha_;

    // effective saturation of phase b
    volScalarField Se_;
        
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
            Se_ == (Sb_+SMALL-Smin_)/(Smax_ - Smin_);
            pc_ = pc0_ * pow(Se_,-alpha_);
            dpcdS_= -alpha_*pc0_*(pow(Se_,-alpha_-1))/(Smax_-Smin_);
        }

    //- Compute Moisture from head pressure (Richards' model)
    volScalarField correctAndSb(volScalarField h) 
        {
            FatalErrorIn("volScalarField MoistureFromHeadPressure(volScalarField h) const")
                << "Not implemented in pcBrooksAndCorey.H " << abort(FatalError);
            return Se_;
        };

    //- Capillary capacity (Richards' model)
    tmp<volScalarField> Ch() const
        {
            FatalErrorIn("tmp<volScalarField> Ch() const")
                << "Not implemented in pcBrooksAndCorey.H " << abort(FatalError);
            return Se_;
        };

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
