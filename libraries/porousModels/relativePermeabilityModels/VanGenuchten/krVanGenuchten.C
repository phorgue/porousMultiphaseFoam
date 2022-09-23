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

#include "krVanGenuchten.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krVanGenuchten, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krVanGenuchten,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krVanGenuchten::krVanGenuchten
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname
)
    :
    relativePermeabilityModel(mesh, transportProperties, Sname),
    krVanGenuchtenCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
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
        dimensionedScalar("m",dimless,krVanGenuchtenCoeffs_.getOrDefault<scalar>("m",0))
    ),
    kramax_
    (
        IOobject
        (
            "kr"+Sname+"max",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kr"+Sname+"max",dimless,krVanGenuchtenCoeffs_.getOrDefault<scalar>("kr"+Sname+"max",1.0))
    ),
    krbmax_
    (
        IOobject
        (
            "kr"+Sname+"max",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kr"+Sname+"max",dimless,krVanGenuchtenCoeffs_.getOrDefault<scalar>("kr"+Sname+"max",1.0))
    )
{
    if (gMin(m_) <= 0)
    {
        FatalErrorIn
            (
                "in krVanGenuchten.C"
            )
            << "Relative permeability coefficient m equal or less than 0" 
                << exit(FatalError);
    }
    dimensionedScalar Smin = krVanGenuchtenCoeffs_.getOrDefault(Sname+"min",dimensionedScalar(Sname+"min", dimless, 0));
    if (Smin.value() > 0) Smin_ = Smin;
    dimensionedScalar Smax = krVanGenuchtenCoeffs_.getOrDefault(Sname+"max",dimensionedScalar(Sname+"max", dimless, 1));
    if (Smax.value() < 1) Smax_ = Smax;
    Info << "Van Genuchten parameters for relative permeability model" << nl << "{" << endl;
    Info << "    m ";
    if (m_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(m_).value() << endl;}
    Info << "    Smax ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info <<  "    Smin ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    kramax ";
    if (kramax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(kramax_).value() << endl;}
    Info << "    krbmax ";
    if (krbmax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(krbmax_).value() << endl;}
    Info << "} \n" << endl;   
}

// ************************************************************************* //
