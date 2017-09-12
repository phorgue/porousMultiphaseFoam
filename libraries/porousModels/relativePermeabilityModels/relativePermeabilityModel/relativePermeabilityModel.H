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
    Foam::relativePermeabilityModel

Description
    A namespace for various relative permeability implementations.

Class
    Foam::relativePermeabilityModel

Description
    An abstract base class for relative permeability models.

SourceFiles
    relativePermeabilityModel.C
    relativePermeabilityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef relativePermeabilityModel_H
#define relativePermeabilityModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class relativePermeabilityModel Declaration
\*---------------------------------------------------------------------------*/

class relativePermeabilityModel
{
protected:

    // Protected data
    word name_;
    dictionary transportProperties_;

    //- reference to saturation field
    const volScalarField& Sb_;

    //- relative permeability and derivative for each phase
    volScalarField kra_,krb_;
    volScalarField dkradS_,dkrbdS_;

    //- Disallow copy construct
    relativePermeabilityModel(const relativePermeabilityModel&);

    //- Disallow default bitwise assignment
    void operator=(const relativePermeabilityModel&);


public:

    //- Runtime type information
    TypeName("relativePermeabilityModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        relativePermeabilityModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const volScalarField& Sb
        ),
        (name, transportProperties, Sb)
    );


    // Selectors

    //- Return a reference to the selected relativePermeability model
    static autoPtr<relativePermeabilityModel> New
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sb
    );


    // Constructors

    //- Construct from components
    relativePermeabilityModel
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sb
    );


    //- Destructor
    virtual ~relativePermeabilityModel()
        {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& transportProperties() const
        {
            return transportProperties_;
        }

    //- Return relative permeabilities according
    tmp<volScalarField> kra() const
        {
            return kra_;
        }

    tmp<volScalarField> krb() const
        {
            return krb_;
        }

    //- return the derivatives of relative permeabilities
    tmp<volScalarField> dkradS() const
        {
            return dkradS_;
        }   
    tmp<volScalarField> dkrbdS() const
        {
            return dkrbdS_;
        }
      
    //- Correct the relative permeabilities
    virtual void correct() = 0; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
