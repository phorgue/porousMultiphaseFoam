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
    Foam::porousMediumModel

Description
    A namespace for various relative permeability implementations.

Class
    Foam::porousMediumModel

Description
    An abstract base class for relative permeability models.

SourceFiles
    porousMediumModel.C
    porousMediumModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef porousMediumModel_H
#define porousMediumModel_H

#include "surfaceFields.H"
#include "volFields.H"
#include "incompressiblePhase.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class porousMediumModel Declaration
\*---------------------------------------------------------------------------*/

class porousMediumModel
{
protected:

    // Protected data
    //- Name of the main saturation field
    const word Sname_;

    //- Reference to the main dictionary
    const dictionary& transportProperties_;

    //- Reference to the phase properties
    const autoPtr<incompressiblePhase>& phase_;

    //- kr/pc models
    autoPtr<relativePermeabilityModel> krModel_;
    autoPtr<capillarityModel> pcModel_;

    //- flow source term
    volScalarField sourceTerm_;

    //- Disallow copy construct
    porousMediumModel(const porousMediumModel&);

    //- Disallow default bitwise assignment
    void operator=(const porousMediumModel&);

public:

    //- Runtime type information
    TypeName("porousMediumModel");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        porousMediumModel,
        dictionary,
        (
            const word Sname,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const autoPtr<incompressiblePhase>& phase,
            const word mediumName
        ),
        (Sname, mesh, transportProperties, phase)
    );

    // Selectors

    //- Return a reference to the selected porousMedium model
    static autoPtr<porousMediumModel> New
    (
        const word Sname,
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const autoPtr<incompressiblePhase>& phase,
        const word mediumName = ""
    );

    // Constructors

    //- Construct from components
    porousMediumModel
    (
        const word Sname,
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const autoPtr<incompressiblePhase>& phase,
        const word mediumName = ""
    );

    //- Destructor
    virtual ~porousMediumModel()
        {}

    // Member Functions
    const word& Sname() const {return Sname_;}

    //- Return transport properties
    const dictionary& transportProperties() const {return transportProperties_;}

    //- Return phase
    const autoPtr<incompressiblePhase>& phase() const {return phase_;}

    //- Return relative permeabilities model
    autoPtr<relativePermeabilityModel>& krModel() {return krModel_;}
    //- Return capillary model
    autoPtr<capillarityModel>& pcModel() {return pcModel_;}

    const volScalarField& source() {return sourceTerm_;}
    
    //- Correct
    virtual void correct() = 0;
    virtual void correct(volScalarField& h, const bool steady = false, const bool conservative = true) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
