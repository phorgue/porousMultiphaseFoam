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

#include "pcIppisch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{
defineTypeNameAndDebug(pcIppisch, 0);

addToRunTimeSelectionTable
(
    capillarityModel,
    pcIppisch,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcIppisch::pcIppisch
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname
)
    :
    capillarityModel(mesh, transportProperties, Sname),
    pcIppischCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    Smin_
    (
        IOobject
        (
            Sname+"min",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        transportProperties.getOrDefault(Sname+"min",dimensionedScalar(Sname+"min",dimless,0))
    ),
    Smax_
    (
        IOobject
        (
            Sname+"max",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        transportProperties.getOrDefault(Sname+"max",dimensionedScalar(Sname+"min",dimless,0))
    ),
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
        dimensionedScalar("m",dimless,pcIppischCoeffs_.getOrDefault<scalar>("m",0))
    ),
    n_(1/(1-m_)),
    alpha_
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
        dimensionedScalar("alpha",dimless,pcIppischCoeffs_.getOrDefault<scalar>("alpha",GREAT))
    ),
    tau_
    (
        IOobject
        (
            "tau",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("tau",dimless,pcIppischCoeffs_.getOrDefault<scalar>("tau",1.))
    ),
    he_
    (
        IOobject
        (
            "he",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("he",dimless,pcIppischCoeffs_.getOrDefault<scalar>("he",1.))
    ),
    Sc_(pow(1+pow(alpha_*he_,n_),-m_))
{
    if (gMin(m_) == 0) FatalErrorIn("Foam::capillarityModels::pcIppisch::pcIppisch") << "m = 0 in pcIppisch" << abort(FatalError);
    Info << "Ippisch parameters for capillary pressure model" << nl << "{" << endl;
    Info << "    m ";
    if (m_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(m_).value() << endl;}
    Info <<  "    n ";
    if (n_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(n_).value() << endl;}
        Info <<  "    alpha ";
    if (alpha_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alpha_).value() << endl;}
        Info <<  "    tau ";
    if (tau_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(tau_).value() << endl;}
    Info <<  "    he ";
    if (he_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(he_).value() << endl;}
    Info << "} \n" << endl;     
}

// ************************************************************************* //

