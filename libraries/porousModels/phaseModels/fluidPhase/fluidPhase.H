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
    Foam::fluidPhase

Description
    Generic fluid phase model for the "porousMultiphaseFoam" toolbox.

SourceFiles
    fluidPhase.C

\*---------------------------------------------------------------------------*/

    #ifndef fluidPhase_H
    #define fluidPhase_H

    #include "phase.H"
    #include "dimensionedScalar.H"
    #include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------* \
                             Class fluidPhase Declaration
    \*---------------------------------------------------------------------------*/

class fluidPhase
:
    public phase
{
   
protected:

    //- phase velocity
    volVectorField U_;

public:

    // Constructors

    fluidPhase
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& phaseName
    );


    // Selectors

    //- Return a reference to the selected turbulence model
    static autoPtr<fluidPhase> New
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& phaseName
    );

    //- Destructor
    virtual ~fluidPhase();

    // Member Functions
    const word& name() const
        {
            return name_;
        }
    
    const volVectorField& U() const
        {
            return U_;
        }

    volVectorField& U()
        {
            return U_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
