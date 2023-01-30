/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

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
    Foam::porousMediumTransportModel

Description
    A namespace for various transport models in porous media implementations.

Class
    Foam::porousMediumTransportModel

Description
    An abstract base class for transport models in porous medium (simple or
    dual porosity)

SourceFiles
    porousMediumTransportModel.C
    porousMediumTransportModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef porousMediumTransportModel_H
#define porousMediumTransportModel_H

#include "surfaceFields.H"
#include "volFields.H"
#include "porousMediumModel.H"
#include "multiscalarMixture.H"
#include "sourceEventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class porousMediumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class porousMediumTransportModel
{
protected:

    // Protected data

    //- Reference to the main dictionary
    const IOdictionary& transportProperties_;

    //- Name of the saturation field for transported species
    word phaseName_;

    //- List of the transported species
    wordList speciesNames_;

    //- Reference to the 
    const porousMediumModel& pmModel_;

    //- List of sourceEvent for each component
    List<sourceEventFile*> sourceEventList_;

    //- Composition
    multiscalarMixture composition_;

    //- Disallow copy construct
    porousMediumTransportModel(const porousMediumTransportModel&);

    //- Disallow default bitwise assignment
    void operator=(const porousMediumTransportModel&);

public:

    //- Runtime type information
    TypeName("porousMediumTransportModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        porousMediumTransportModel,
        dictionary,
        (
            const word& phaseName,
            const porousMediumModel& pmModel
        ),
        (phaseName, pmModel)
    );

    // Selectors

    //- Return a reference to the selected porousMedium model
    static autoPtr<porousMediumTransportModel> New
    (
        const word& phaseName,
        const porousMediumModel& pmModel
    );

    // Constructors

    //- Construct from components
    porousMediumTransportModel
    (
        const word& phaseName,
        const porousMediumModel& pmModel
    );

    //- Destructor
    virtual ~porousMediumTransportModel()
        {}

    // Member Functions

    //- Return transport properties
    const dictionary& transportProperties() const {return transportProperties_;}

    //- Return species names (adding option porous region name)
    wordList speciesNames(word porousRegion="");

    //- Return porousMediumModel
    const porousMediumModel& pmModel() const {return pmModel_;}

    //- Access to the list of sourceEvent
    List<sourceEventFile*>& sourceEventList() {return sourceEventList_;}

    //- Access to composition
    multiscalarMixture& composition() {return composition_;}

    //- Solve tracer tranport equation
    virtual void solveTransport
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& theta
    ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
