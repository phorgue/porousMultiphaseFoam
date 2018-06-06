/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::fixedHeadPressure

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides head pressure condition for h,
    calculated as:

        \f[
            h = potential_ - z (if potential > z)
            h = hUnstaturated (if potential < z)
        \f]

    where
    \vartable
        h       | Head pressure [m]
        hFace   | Height in the opposite direction to gravity
        hRef    | Reference height in the opposite direction to gravity
        \rho    | density
        g       | acceleration due to gravity [m/s^2]
    \endtable

Usage
    \table
        Property     | Description                  | Required    | Default value
        potential    | height of h=0                | no          | 0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedHeadPressure;
        potential       50;
        hUnsaturated    -0.87;
        value           uniform 0; // optional initial value
    }
    \endverbatim

See also
    Foam::fixedValue

SourceFiles
    fixedHeadPressure.C

\*---------------------------------------------------------------------------*/

#ifndef fixedHeadPressure_H
#define fixedHeadPressure_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class fixedHeadPressure Declaration
\*---------------------------------------------------------------------------*/

class fixedHeadPressure
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data
    
    //- Fixed potential
    scalar potential_;

public:

    //- Runtime type information
    TypeName("fixedHeadPressure");


    // Constructors

        //- Construct from patch and internal field
        fixedHeadPressure
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedHeadPressure
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedHeadPressure onto a new patch
        fixedHeadPressure
        (
            const fixedHeadPressure&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedHeadPressure
        (
            const fixedHeadPressure&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new fixedHeadPressure(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedHeadPressure
        (
            const fixedHeadPressure&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedHeadPressure(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
