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
    const word& Sname,
    const word mediumName
)
    :
    capillarityModel(mesh, transportProperties, Sname, mediumName),
    pcIppischCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    m_
    (
        IOobject
        (
            "m"+mediumName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,pcIppischCoeffs_.getOrDefault<scalar>("m"+mediumName,0))
    ),
    n_(1/(1-m_)),
    alpha_
    (
        IOobject
        (
            "alpha"+mediumName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,pcIppischCoeffs_.getOrDefault<scalar>("alpha"+mediumName,GREAT))
    ),
    tau_
    (
        IOobject
        (
            "tau"+mediumName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,pcIppischCoeffs_.getOrDefault<scalar>("tau"+mediumName,1.))
    ),
    he_
    (
        IOobject
        (
            "he"+mediumName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,pcIppischCoeffs_.getOrDefault<scalar>("he"+mediumName,1.))
    ),
    Sc_(pow(1+pow(alpha_*he_,n_),-m_))
{
    dimensionedScalar Smin = pcIppischCoeffs_.getOrDefault(Sname+"min",dimensionedScalar(Sname+"min", dimless, 0));
    if (Smin.value() > 0) Smin_ = Smin;
    dimensionedScalar Smax = pcIppischCoeffs_.getOrDefault(Sname+"max",dimensionedScalar(Sname+"max", dimless, 1));
    if (Smax.value() < 1) Smax_ = Smax;
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

