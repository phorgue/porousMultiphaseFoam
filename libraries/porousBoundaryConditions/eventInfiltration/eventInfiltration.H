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
    Foam::eventInfiltration

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a time varying fixed value condition 
    for velocity field computed using a constant infiltration value 
    given as parameter and variable user-defined patch event file.

Usage
    \table
        Property              | Description                  | Required    | Default value
        constantValue         | constant infiltration value  | no          | 0
        eventFile             | File name                    | no          | ""
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                  eventInfiltration;
        constantValue         -0.01;
        eventFile             "patchEventFileExample.dat";
        value                 uniform 0; // value for paraview visualization
    }
    \endverbatim

See also
    Foam::fixedValue

SourceFiles
    eventInfiltration.C

\*---------------------------------------------------------------------------*/

#ifndef eventInfiltration_H
#define eventInfiltration_H

#include "fixedValueFvPatchFields.H"
#include "patchEventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class eventInfiltration Declaration
\*---------------------------------------------------------------------------*/

class eventInfiltration
:
    public fixedValueFvPatchVectorField
{

protected:

    // Protected data
    
    //- Fixed flux value
    scalar fixedInfiltrationValue_;

    //- time scheme selector
    bool isBackwardScheme_;

    //- Event list (optional)
    label patchEventID_;
    patchEventFile eventFile_;

public:

    //- Runtime type information
    TypeName("eventInfiltration");


    // Constructors

        //- Construct from patch and internal field
        eventInfiltration
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        eventInfiltration
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  eventInfiltration onto a new patch
        eventInfiltration
        (
            const eventInfiltration&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        eventInfiltration
        (
            const eventInfiltration&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField >
            (
                new eventInfiltration(*this)
            );
        }

        //- Construct as copy setting internal field reference
        eventInfiltration
        (
            const eventInfiltration&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new eventInfiltration(*this, iF)
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
