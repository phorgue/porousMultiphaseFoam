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
    Foam::phase

Description
    Generic phase model for the "porousMultiphaseFoam" toolbox.

SourceFiles
    phase.C

\*---------------------------------------------------------------------------*/

    #ifndef phase_H
    #define phase_H

    #include "dictionary.H"
    #include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------* \
                             Class fluidPhase Declaration
    \*---------------------------------------------------------------------------*/

class phase
{

protected:

    //- reference to the mesh
    const fvMesh& mesh_;

    // Private data
    dictionary dict_;

    //- Name of phase
    word name_;

public:

    // Constructors

    phase
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& phaseName
    );

    // Selectors

    //- Return a reference to the selected turbulence model
    static autoPtr<phase> New
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& phaseName
    );

    //- Destructor
    virtual ~phase();

    // Member Functions
    const word& name() const
        {
            return name_;
        }

    const fvMesh& mesh() const
        {
            return mesh_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
