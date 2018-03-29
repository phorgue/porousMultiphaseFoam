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
    Foam::fixedHeadPressureSTL

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides head pressure condition for h,
    calculated as the distance between the face and the closest point
    of the given surface :

        \f[
            h = zSurface - zFace
        \f]

    where
    \vartable
        h       | Head pressure [m]
        zSurface| Local height of the STL surface
        zFace   | Height of the face center
    \endtable

Usage
    \table
        Property     | Description                  | Required    | Default value
        STL          | pattern to the surface file  | yes          | 0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedHeadPressureSTL;
        potential       constant/triSurface/example.stl;
        value           uniform 0; // optional initial value
    }
    \endverbatim

See also
    Foam::fixedValue

SourceFiles
    fixedHeadPressureSTL.C

\*---------------------------------------------------------------------------*/

#ifndef fixedHeadPressureSTL_H
#define fixedHeadPressureSTL_H

#include "fixedValueFvPatchFields.H"
 #include "triSurfaceMesh.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class fixedHeadPressureSTL Declaration
\*---------------------------------------------------------------------------*/

class fixedHeadPressureSTL
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data

    //- Name of phase-fraction field
    word STLname_;
    
public:

    //- Runtime type information
    TypeName("fixedHeadPressureSTL");


    // Constructors

        //- Construct from patch and internal field
        fixedHeadPressureSTL
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedHeadPressureSTL
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedHeadPressureSTL onto a new patch
        fixedHeadPressureSTL
        (
            const fixedHeadPressureSTL&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedHeadPressureSTL
        (
            const fixedHeadPressureSTL&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new fixedHeadPressureSTL(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedHeadPressureSTL
        (
            const fixedHeadPressureSTL&,
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
                new fixedHeadPressureSTL(*this, iF)
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
