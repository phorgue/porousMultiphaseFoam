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

Class
    Foam::relativePermeabilityModels::krVanGenuchten

Description
    Standard Van Genuchten relative permeability model.

SourceFiles
    krVanGenuchten.C

\*---------------------------------------------------------------------------*/

#ifndef krVanGenuchten_H
#define krVanGenuchten_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
    Class krVanGenuchten Declaration
    \*---------------------------------------------------------------------------*/

class krVanGenuchten
:
    public relativePermeabilityModel
{
    //- Van Genuchten coefficients
    volScalarField m_;

    //- end points
    volScalarField kramax_;
    volScalarField krbmax_;

public:

    //- Runtime type information
    TypeName("VanGenuchten");

    // Constructors

    //- Construct from components
    krVanGenuchten
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~krVanGenuchten()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct(const volScalarField& Sb, bool derivative);
    void correctkra(const volScalarField& Sb, bool derivative);
    void correctkrb(const volScalarField& Sb, bool derivative);
    void correctkrb(const volScalarField& Sb, const label& celli);
    tmp<volScalarField> kr(const volScalarField& S);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
