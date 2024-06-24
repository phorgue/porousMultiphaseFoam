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

#include "dualAlphaDispersion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
defineTypeNameAndDebug(dualAlphaDispersion, 0);

addToRunTimeSelectionTable
(
    dispersionModel,
    dualAlphaDispersion,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModels::dualAlphaDispersion::dualAlphaDispersion
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
    :
    dispersionModel(name, transportProperties, mesh),
    dualAlphaDispersionCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    tau_
    (
        IOobject
        (
            "tau",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("tau",dualAlphaDispersionCoeffs_)
    ),
    alphaLsat_
    (
        IOobject
        (
            "alphaL_sat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaL_sat",dualAlphaDispersionCoeffs_)
    ),
    alphaTsat_
    (
        IOobject
        (
            "alphaT_sat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaT_sat",dualAlphaDispersionCoeffs_)
    ),
    alphaLunsat_
    (
        IOobject
        (
            "alphaL_unsat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaL_unsat",dualAlphaDispersionCoeffs_)
    ),
    alphaTunsat_
    (
        IOobject
        (
            "alphaT_unsat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaT_unsat",dualAlphaDispersionCoeffs_)
    ),
    coef_
    (
        IOobject
        (
            "coef",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("coef", dimless, 0)
    ),
    thetaThresholdSat_
    (
        IOobject
        (
            "thetaThreshold_sat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("thetaThreshold_sat",dualAlphaDispersionCoeffs_)
    ),
    thetaThresholdUnsat_
    (
        IOobject
        (
            "thetaThreshold_unsat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("thetaThreshold_unsat",dualAlphaDispersionCoeffs_)
    )  
{
    Info << "Parameters for dualAlpha dispersion model" << nl << "{" << endl;
    Info << "    tau ";
    if (tau_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(tau_).value() << endl;}
    Info << "    alphaL_sat ";
    if (alphaLsat_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alphaLsat_).value() << endl;}
    Info << "    alphaT_sat ";
    if (alphaTsat_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alphaTsat_).value() << endl;}
    Info <<  "   alphaL_unsat ";
    if (alphaLunsat_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alphaLunsat_).value() << endl;}
    Info << "    alphaT_unsat ";
    if (alphaTunsat_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alphaTunsat_).value() << endl;}
    Info << "    thetaThreshold_sat ";
    if (thetaThresholdSat_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(thetaThresholdSat_).value() << endl;}
    Info << "    thetaThreshold_unsat ";
    if (thetaThresholdUnsat_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(thetaThresholdUnsat_).value() << endl;}
    Info << "}" << endl;
}

// ************************************************************************* //
