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
    Foam::porousMediumModels::simplePorosity

Class
    Foam::simplePorosity

Description
    Simple porosity model (one kr / one pc law)

SourceFiles
    simplePorosity.C

\*---------------------------------------------------------------------------*/

#ifndef simplePorosity_H
#define simplePorosity_H

#include "volFields.H"
#include "porousMediumModel.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumModels
{

/*---------------------------------------------------------------------------* \
             Class simplePorosity Declaration
\*---------------------------------------------------------------------------*/

class simplePorosity
:
    public porousMediumModel
{
protected:

    //- Disallow copy construct
    simplePorosity(const simplePorosity&);

    //- Disallow default bitwise assignment
    void operator=(const simplePorosity&);

public:

    //- Runtime type information
    TypeName("simplePorosity");

    // Constructors

    //- Construct from components
    simplePorosity
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        autoPtr<relativePermeabilityModel>& krModel,
        autoPtr<capillarityModel>& pcModel
    );

    //- Destructor
    ~simplePorosity()
        {}

    //- Correct the relative permeabilities
    void correct();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousMediumModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
