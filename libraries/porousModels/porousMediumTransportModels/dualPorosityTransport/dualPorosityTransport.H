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
    Foam::porousMediumModels::dualPorosityTransport

Class
    Foam::dualPorosityTransport

Description
    Dual porosity model (one kr / one pc law)

SourceFiles
    dualPorosityTransport.C

\*---------------------------------------------------------------------------*/

#ifndef dualPorosityTransport_H
#define dualPorosityTransport_H

#include "volFields.H"
#include "surfaceFields.H"
#include "porousMediumTransportModel.H"
#include "incompressiblePhase.H"
#include "multiscalarMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumTransportModels
{

/*---------------------------------------------------------------------------* \
             Class dualPorosityTransport Declaration
\*---------------------------------------------------------------------------*/

class dualPorosityTransport
:
    public porousMediumTransportModel
{
    //- Dictionary for matrix parameters
    const dictionary dualPorosityTransportCoeffs_;

    //- List of sourceEvent for each component
    List<sourceEventFile*> sourceEventList_;

    //- Porosity for matrix retard coefficient
    volScalarField epsMatrix_;

    //- Composition
    multiscalarMixture matrixComposition_;

    //- Fracture/matrix echange coefficent
    dimensionedScalar a_, beta_, gammaW_, alphaS_;

    //- Fracture/matrix fields for tracer-water transfers
    volScalarField exchangeTermFromFracture_, exchangeTermFromMatrix_;

    //- Reference to required fields
    const volVectorField& UMatrix_;
    const surfaceScalarField& phiMatrix_;
    const volScalarField& thetaMatrix_;

    //- Disallow copy construct
    dualPorosityTransport(const dualPorosityTransport&);

    //- Disallow default bitwise assignment
    void operator=(const dualPorosityTransport&);

public:

    //- Runtime type information
    TypeName("dualPorosityTransport");

    // Constructors

    //- Construct from components
    dualPorosityTransport
    (
        const word& phaseName,
        const porousMediumModel& pmModel
    );

    //- Destructor
    ~dualPorosityTransport()
        {}

    //- Solve tracer tranport equation
    void solveTransport
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& theta
    );

    //- Display mass balance information to the terminal
    void info
    (
        const label& speciei,
        const volScalarField& theta,
        const surfaceScalarField& phi
    ) ;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousMediumTransportModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
