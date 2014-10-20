/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::incompressiblePhase

Description
    Incompressible phase model for the "porousMultiphaseFoam" toolbox.

SourceFiles
    incompressiblePhase.C

\*---------------------------------------------------------------------------*/

    #ifndef incompressiblePhase_H
    #define incompressiblePhase_H

    #include "dictionary.H"
    #include "dimensionedScalar.H"
    #include "volFields.H"
    #include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                   Class incompressiblePhase Declaration
    \*---------------------------------------------------------------------------*/

    class incompressiblePhase
    {
        //- Name of phase
        word name_;

        // Private data
        dictionary dict_;

        //- dynamic viscosity
        dimensionedScalar mu_;

        //- density
        dimensionedScalar rho_;

        //- phase velocity
        volVectorField U_;

        //- phase fluxes
        autoPtr<surfaceScalarField> phiPtr_;

    public:

        // Constructors

        incompressiblePhase
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& phaseName
        );


        // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<incompressiblePhase> New
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& phaseName
        );

        //- Destructor
        virtual ~incompressiblePhase();

        // Member Functions
        const word& name() const
            {
                return name_;
            }

        const dimensionedScalar& mu() const
            {
                return mu_;
            }

        const dimensionedScalar& rho() const
            {
                return rho_;
            }
      
        const volVectorField& U() const
            {
                return U_;
            }

        volVectorField& U()
            {
                return U_;
            }

        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        surfaceScalarField& phi()
        {
            return phiPtr_();
        }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
