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

\*---------------------------------------------------------------------------*/

#include "pcVanGenuchten.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{
defineTypeNameAndDebug(pcVanGenuchten, 0);

addToRunTimeSelectionTable
(
    capillarityModel,
    pcVanGenuchten,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcVanGenuchten::pcVanGenuchten
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname
)
    :
    capillarityModel(mesh, transportProperties, Sname),
    pcVanGenuchtenCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    m_
    (
        IOobject
        (
            "m",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("m",dimless,pcVanGenuchtenCoeffs_.getOrDefault<scalar>("m",0))
    ),
    n_(1/(1-m_)),
    alpha_ // necessary for Richards solver
    (
        IOobject
        (
            "alpha",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha",dimless,pcVanGenuchtenCoeffs_.getOrDefault<scalar>("alpha",0))
    ),
    pc0_
    (
        IOobject
        (
            "pc0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        pcVanGenuchtenCoeffs_.getOrDefault("pc0",dimensionedScalar("pc0",dimensionSet(1,-1,-2,0,0),0.))
    )
{
    dimensionedScalar Smin = pcVanGenuchtenCoeffs_.getOrDefault(Sname+"min",dimensionedScalar(Sname+"min", dimless, 0));
    if (Smin.value() > 0) Smin_ = Smin;
    dimensionedScalar Smax = pcVanGenuchtenCoeffs_.getOrDefault(Sname+"max",dimensionedScalar(Sname+"max", dimless, 1));
    if (Smax.value() < 1) Smax_ = Smax;
    if (gMin(m_) == 0) FatalErrorIn("Foam::capillarityModels::pcVanGenuchten::pcVanGenuchten") << "m = 0 in pcVanGenuchten" << abort(FatalError);
    Info << "Van Genuchten parameters for capillary pressure model" << nl << "{" << endl;
    Info << "    m ";
    if (m_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(m_).value() << endl;}
    Info << "    pc0 ";
    if (pc0_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(pc0_).value() << endl;}
    Info << "    alpha ";
    if (alpha_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alpha_).value() << endl;}
    Info <<  "    Smin ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    Smax ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info << "} \n" << endl;
    
}

// ************************************************************************* //
