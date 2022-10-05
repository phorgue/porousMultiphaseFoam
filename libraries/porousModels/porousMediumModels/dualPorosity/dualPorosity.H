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
    Foam::porousMediumModels::dualPorosity

Class
    Foam::dualPorosity

Description
    Dual porosity model (one kr / one pc law)

SourceFiles
    dualPorosity.C

\*---------------------------------------------------------------------------*/

#ifndef dualPorosity_H
#define dualPorosity_H

#include "volFields.H"
#include "surfaceFields.H"
#include "porousMediumModel.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"
#include "incompressiblePhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumModels
{

/*---------------------------------------------------------------------------* \
             Class dualPorosity Declaration
\*---------------------------------------------------------------------------*/

class dualPorosity
:
    public porousMediumModel
{
    const word SnameMatrix_;
    const fvMesh& mesh_;
    const uniformDimensionedVectorField g;
    volScalarField hMatrix_;
    volScalarField Smatrix_;
    volScalarField KMatrix_;
    surfaceScalarField KMatrixf_;
    volScalarField Wf_;
    volScalarField  geomFactor_;
    volVectorField UMatrix_;
    surfaceScalarField LMatrixf_;
    surfaceScalarField MMatrixf_;
    surfaceScalarField phiGMatrixf_;
    surfaceScalarField phiMatrix_;

    //- kr/pc model for the matrix
    autoPtr<relativePermeabilityModel> matrixKrModel_;
    autoPtr<capillarityModel> matrixPcModel_;

    //- Update properties of the matrix part
    void updateMatrixProperties();

    //- Disallow copy construct
    dualPorosity(const dualPorosity&);

    //- Disallow default bitwise assignment
    void operator=(const dualPorosity&);

public:

    //- Runtime type information
    TypeName("dualPorosity");

    // Constructors

    //- Construct from components
    dualPorosity
    (
        const word Sname,
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const autoPtr<incompressiblePhase>& phase
    );

    //- Destructor
    ~dualPorosity()
        {}

    //- Correct the relative permeabilities
    void correct();
    void correct(volScalarField& h, const bool steady = false, const bool conservative = true);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousMediumModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
