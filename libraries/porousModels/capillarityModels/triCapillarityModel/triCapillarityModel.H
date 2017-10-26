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
    Foam::triCapillarityModel

Description
    A namespace for various 3-phases capillarity models.

Class
    Foam::triCapillarityModel

Description
    An abstract base class for 3-phases capillarity models.

SourceFiles
    triCapillarityModel.C
    triCapillarityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef triCapillarityModel_H
#define triCapillarityModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class triCapillarityModel Declaration
\*---------------------------------------------------------------------------*/

class triCapillarityModel
{

protected:

    // Protected data

    word name_;
    dictionary triCapillarityProperties_;

    // reference to the saturation field 
    const volScalarField& Sw_;
    const volScalarField& So_;
    
    // capillary pressure and derivative
    volScalarField pcao_,pcow_,pcaw_,dpcaodS_,dpcowdS_,dpcawdS_;

    //- Disallow copy construct
    triCapillarityModel(const triCapillarityModel&);

    //- Disallow default bitwise assignment
    void operator=(const triCapillarityModel&);

public:

    //- Runtime type information
    TypeName("triCapillarityModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        triCapillarityModel,
        dictionary,
        (
            const word& name,
            const dictionary& triCapillarityProperties,
            const volScalarField& Sw,
            const volScalarField& So
        ),
        (name, triCapillarityProperties, Sw, So)
    );


    // Selectors

    //- Return a reference to the selected triCapillarity model
    static autoPtr<triCapillarityModel> New
    (
        const word& name,
        const dictionary& triCapillarityProperties,
        const volScalarField& Sw,
        const volScalarField& So
    );


    // Constructors

    //- Construct from components
    triCapillarityModel
    (
        const word& name,
        const dictionary& triCapillarityProperties,
        const volScalarField& Sw,
        const volScalarField& So
    );


    //- Destructor
    virtual ~triCapillarityModel() {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& triCapillarityProperties() const
        {
            return triCapillarityProperties_;
        }

    //- Return the capillary pressure according to Parker and Lenhard correlation
    tmp<volScalarField> pcao() const
        {
            return pcao_;
        }
        
    tmp<volScalarField> pcow() const
        {
            return pcow_;
        }

    tmp<volScalarField> pcaw() const
        {
            return pcaw_;
        }
    
    //- Return the capillary pressure derivative according to Parker and Lenhard correlation
    tmp<volScalarField> dpcaodS() const
        {
            return dpcaodS_;
        }
        
    tmp<volScalarField> dpcowdS() const
        {
            return dpcowdS_;
        }
    
    tmp<volScalarField> dpcawdS() const
        {
            return dpcawdS_;
        }
    
    //- Correct the 3-phases capillary pressures
    virtual void correct() = 0; 

    //- Compute Moisture from head pressure (Richards' model)
    //virtual volScalarField correctAndSb(volScalarField h) = 0;

    //- TriCapillary capacity (Richards' model)
    //virtual tmp<volScalarField> Ch() const = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
