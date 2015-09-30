/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::capillarityModel

Description
    A namespace for various capillarity Effect implementations.

Class
    Foam::capillarityModel

Description
    An abstract base class for incompressible capillarityModels.

SourceFiles
    capillarityModel.C
    capillarityModelNew.C

    \*---------------------------------------------------------------------------*/

#ifndef capillarityModel_H
#define capillarityModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class capillarityModel Declaration
    \*---------------------------------------------------------------------------*/

class capillarityModel
{

protected:

    // Protected data

    word name_;
    dictionary capillarityProperties_;

    // reference to the saturation field 
    const volScalarField& Sb_;

    // Private Member Functions

    //- Disallow copy construct
    capillarityModel(const capillarityModel&);

    //- Disallow default bitwise assignment
    void operator=(const capillarityModel&);

public:

    //- Runtime type information
    TypeName("capillarityModel");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        capillarityModel,
        dictionary,
        (
            const word& name,
            const dictionary& capillarityProperties,
            const volScalarField& Sb
        ),
        (name, capillarityProperties, Sb)
    );


    // Selectors

    //- Return a reference to the selected capillarity model
    static autoPtr<capillarityModel> New
    (
        const word& name,
        const dictionary& capillarityProperties,
        const volScalarField& Sb
    );


    // Constructors

    //- Construct from components
    capillarityModel
    (
        const word& name,
        const dictionary& capillarityProperties,
        const volScalarField& Sb
    );


    //- Destructor
    virtual ~capillarityModel() {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& capillarityProperties() const
    {
      return capillarityProperties_;
    }
      
    //- Return the capillary pressure as a function of Sb 
    virtual tmp<volScalarField> pc() const = 0;

    //- Return the derivative of capillary pressure as a function of Sb
    virtual tmp<volScalarField> dpcdS() const = 0;
        
    //- Correct the capillary pressure
    virtual void correct() = 0; 

    //- Compute Moisture from head pressure (Richards' model)
    virtual volScalarField correctAndSb(volScalarField h) = 0;

    //- Capillary capacity (Richards' model)
    virtual tmp<volScalarField> Ch() const = 0;

    //- Read transportProperties dictionary
    virtual bool read(const dictionary& capillarityProperties) = 0;

};

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
