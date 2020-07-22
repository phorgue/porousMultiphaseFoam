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
    Foam::dispersionModels::alphaDispersion

Description
    Dispersion model based on alphaL and alphaT coefficients (Bear, 1979)

      Deff= Dm/tau + (alphaT * mag(U) + (alphaL - alphaT) * (U*U)/mag(U) )/theta

    where theta is water content :

      theta = Saturation * porosity

SourceFiles
    alphaDispersion.C

Developer
    Pierre Horgue

\*---------------------------------------------------------------------------*/

#ifndef alphaDispersion_H
#define alphaDispersion_H

#include "dispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{

/*---------------------------------------------------------------------------*\
                  Class alphaDispersion Declaration
\*---------------------------------------------------------------------------*/

class alphaDispersion
:
    public dispersionModel
{
    //- Private data

    //- Dictionary for model coefficients
    dictionary alphaDispersionCoeffs_;

    //- model coefficients
    volScalarField tau_,alphaL_,alphaT_;

public:

    //- Runtime type information
    TypeName("alphaDispersion");

    // Constructors

    //- Construct from components
    alphaDispersion
    (
        const word& name,
        const dictionary& transportProperties,
        const fvMesh& mesh
    );

    //- Destructor
    ~alphaDispersion()
        {}

    // Member Functions

    //- Correct the dispersion giving water content or porosity field
    void correct(const volScalarField& C, const volVectorField& U, const volScalarField& theta)
        {
            dimensionedScalar USMALL("",dimLength/dimTime,1e-36);

            Deff_ = (Dm_/tau_ + alphaT_ * mag(U)/theta)*tensor::I + (alphaL_-alphaT_) * (U*U)/(theta*mag(U)+USMALL);
            forAll(Deff_.boundaryField(),patchi)
            {
                //- Dispersion equal to zero on patch where flux is fixed
                if (C.boundaryField().types()[patchi] == "eventFlux" ||
                    C.boundaryField().types()[patchi] == "inletOutlet"
                )
                {
                    Deff_.boundaryFieldRef()[patchi] = tensor(zero());
                }
            }
        }

    //- Correct the dispersion giving scalar porosity
    void correct(const volScalarField& C, const volVectorField& U, const dimensionedScalar& eps)
        {
            dimensionedScalar USMALL("",dimLength/dimTime,1e-36);

            Deff_ = (Dm_/tau_ + alphaT_ * mag(U)/eps)*tensor::I + (alphaL_-alphaT_) * (U*U)/(eps*mag(U)+USMALL);
            forAll(Deff_.boundaryField(),patchi)
            {
                //- Dispersion equal to zero on patches
                if (C.boundaryField().types()[patchi] == "eventFlux" ||
                C.boundaryField().types()[patchi] == "inletOutlet")
                {
                    Deff_.boundaryFieldRef()[patchi] = tensor(zero());
                }
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dispersionModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
