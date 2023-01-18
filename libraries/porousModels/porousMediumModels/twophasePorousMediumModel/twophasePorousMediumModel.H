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
    Foam::twophasePorousMediumModel

Description
    A namespace for various relative permeability implementations.

Class
    Foam::twophasePorousMediumModel

Description
    An abstract base class for relative permeability models.

SourceFiles
    twophasePorousMediumModel.C
    twophasePorousMediumModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef twophasePorousMediumModel_H
#define twophasePorousMediumModel_H

#include "surfaceFields.H"
#include "volFields.H"
#include "porousMediumModel.H"
#include "incompressiblePhase.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class twophasePorousMediumModel Declaration
\*---------------------------------------------------------------------------*/

class twophasePorousMediumModel
:
    public porousMediumModel
{
protected:

    // Protected data
    //- Name of the main saturation field
    const word Sname_;

    //- flow source term
    volScalarField sourceTerm_;

    //- kr/pc models
    autoPtr<relativePermeabilityModel> krModel_;
    autoPtr<capillarityModel> pcModel_;

    //- Disallow copy construct
    twophasePorousMediumModel(const twophasePorousMediumModel&);

    //- Disallow default bitwise assignment
    void operator=(const twophasePorousMediumModel&);

public:

    //- Runtime type information
    TypeName("twophasePorousMediumModel");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        twophasePorousMediumModel,
        dictionary,
        (
            const word Sname,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const autoPtr<incompressiblePhase>& phase,
            const word porousRegion
        ),
        (Sname, mesh, transportProperties, phase, porousRegion)
    );

    // Selectors

    //- Return a reference to the selected twophasePorousMedium model
    static autoPtr<twophasePorousMediumModel> New
    (
        const word Sname,
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const autoPtr<incompressiblePhase>& phase,
        const word porousRegion = ""
    );

    // Constructors

    //- Construct from components
    twophasePorousMediumModel
    (
        const word Sname,
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const autoPtr<incompressiblePhase>& phase,
        const word porousRegion = ""
    );

    //- Destructor
    virtual ~twophasePorousMediumModel()
        {}

    // Member Functions
    const word& Sname() const {return Sname_;}

    //- Return relative permeabilities model
    autoPtr<relativePermeabilityModel>& krModel() {return krModel_;}
    //- Return capillary model
    autoPtr<capillarityModel>& pcModel() {return pcModel_;}

    const volScalarField& source() {return sourceTerm_;}

    //- Correct
    virtual void rewindTime() = 0;
    virtual void correct() = 0;
    virtual void correct(volScalarField& h, const bool steady = false, const bool conservative = true) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
