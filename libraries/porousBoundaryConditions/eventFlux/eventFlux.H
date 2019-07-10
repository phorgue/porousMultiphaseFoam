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
    Foam::eventFlux

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a fixed value condition for a given
    variable C calculated as:

        \f[
            C = (constant value + event value) / totalFlux
        \f]

    where
    \vartable
        C                | Variable solved
        constantValue    | fixed flux value of variable solved
        totalFlux        | computed flux using field name
    \endtable

Usage
    \table
        Property          | Description                  | Required    | Default value
        constantValue     | constant fixed flux          | yes         | 0
        phiName           | Name of the flux field       | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            eventFlux;
        eventFluxValue  50;
        phiName         phi; // optional (phi by default)
        value           uniform 0; // optional initial value
    }
    \endverbatim

See also
    Foam::fixedValue

SourceFiles
    eventFlux.C

\*---------------------------------------------------------------------------*/

#ifndef eventFlux_H
#define eventFlux_H

#include "fixedValueFvPatchFields.H"
#include "patchEventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class eventFlux Declaration
\*---------------------------------------------------------------------------*/

class eventFlux
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data
    
    //- Fixed flux value
    scalar eventFluxValue_;

    //- Flux field name
    word phiName_;
    bool isBackwardScheme_;

    //- Event list (optional)
    label patchEventID_;
    patchEventFile eventFile_;

public:

    //- Runtime type information
    TypeName("eventFlux");


    // Constructors

        //- Construct from patch and internal field
        eventFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        eventFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  eventFlux onto a new patch
        eventFlux
        (
            const eventFlux&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        eventFlux
        (
            const eventFlux&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new eventFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        eventFlux
        (
            const eventFlux&,
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
                new eventFlux(*this, iF)
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
