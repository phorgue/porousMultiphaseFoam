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

Namespace
    Foam::triRelativePermeabilityModel

Description
    A namespace for various 3-phases relative permeability implementations.

Class
    Foam::triRelativePermeabilityModel

Description
    An abstract base class for incompressible 3-phases relative permeability models.

SourceFiles
    triRelativePermeabilityModel.C
    triRelativePermeabilityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef triRelativePermeabilityModel_H
#define triRelativePermeabilityModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class triRelativePermeabilityModel Declaration
    \*---------------------------------------------------------------------------*/

class triRelativePermeabilityModel
{

protected:

    // Protected data
    word name_;
    dictionary transportProperties_;

    //- Reference to saturation fields (water and oil)
    const volScalarField& Sw_;
    const volScalarField& So_;

    //- Relative permeability for each phase
    volScalarField kra_,kro_,krw_;
    volScalarField dkradS_, dkrwdS_;

    // Private Member Functions

    //- Disallow copy construct
    triRelativePermeabilityModel(const triRelativePermeabilityModel&);

    //- Disallow default bitwise assignment
    void operator=(const triRelativePermeabilityModel&);


public:

    //- Runtime type information
    TypeName("triRelativePermeabilityModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        triRelativePermeabilityModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const volScalarField& Sw,
            const volScalarField& So
        ),
        (name, transportProperties, Sw, So)
    );


    // Selectors

    //- Return a reference to the selected triRelativePermeability model
    static autoPtr<triRelativePermeabilityModel> New
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sw,
        const volScalarField& So
    );


    // Constructors

    //- Construct from components
    triRelativePermeabilityModel
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sw,
        const volScalarField& So
    );


    //- Destructor
    virtual ~triRelativePermeabilityModel()
        {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& transportProperties() const
        {
            return transportProperties_;
        }
        
    //- Return the relative permeability of phase a
    tmp<volScalarField> kra() const
        {
            return kra_;
        }

    //- Return the relative permeability of phase c
    tmp<volScalarField> kro() const
        {
            return kro_;
        }
    
    //- Return the relative permeability of phase b
    tmp<volScalarField> krw() const
        {
            return krw_;
        }
    
    //- return the derivatives of relative permeabilities
    tmp<volScalarField> dkradS() const
        {
            return dkradS_;
        }   
    tmp<volScalarField> dkrwdS() const
        {
            return dkrwdS_;
        }
    
    //- Correct the relative permeabilities
    virtual void correct() = 0; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
