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

#include "fvCFD.H"
#include "dualPorosity.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumModels
{
defineTypeNameAndDebug(dualPorosity, 0);

addToRunTimeSelectionTable
(
    porousMediumModel,
    dualPorosity,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumModels::dualPorosity::dualPorosity
(
    const word Sname,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const autoPtr<incompressiblePhase>& phase
)
    :
    porousMediumModel(Sname, mesh, transportProperties, phase),
    SnameM_(Sname+"M"),
    hMatrix_
    (
        IOobject
        (
            "hM",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless,
        calculatedFvPatchScalarField::typeName
    ),
    Smatrix_
    (
        IOobject
        (
            SnameM_,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless,
        calculatedFvPatchScalarField::typeName
    ),
    Km_
    (
        IOobject
        (
            "K",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Kmf_(fvc::interpolate(Km_, "K")),
    Umatrix_
    (
        IOobject
        (
            "UM",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    phiMatrix_
    (
        IOobject
        (
            "phiM",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Umatrix_) & mesh.Sf()
    )
{
    matrixPcModel_ = capillarityModel::New(mesh, transportProperties, SnameM_);
    matrixKrModel_ = relativePermeabilityModel::New(mesh, transportProperties, SnameM_);
}

void Foam::porousMediumModels::dualPorosity::correct()
{
    FatalErrorIn("dualPorosity.C") << " dualPorosity cannot be used with impesFoam/anisoImpesFoam " << abort(FatalError);
}

void Foam::porousMediumModels::dualPorosity::correct(const volScalarField& h, const bool steady, const bool conservative)
{

    const meshObjects::gravity& g = meshObjects::gravity::New(hMatrix_.time());
    surfaceScalarField krthetaMf("krthetaf",fvc::interpolate(matrixKrModel_->krb(),"krtheta"));

    surfaceScalarField Lf ("Lf",phase_->rho()*Kmf_*krthetaMf/phase_->mu());
    surfaceScalarField Mf ("Mf",mag(g)*Lf);
    surfaceScalarField phiG("phiG",(Lf * g) & hMatrix_.mesh().Sf());

    hMatrix_.storePrevIter();
    
    //- solve matrix equation
    fvScalarMatrix hMEqn
        (
            //- transport terms
            - fvm::laplacian(Mf,hMatrix_)
            + fvc::div(phiG)
        );
    if (!steady)
    {
        //- accumulation terms
        hMEqn += matrixPcModel_->Ch() * fvm::ddt(hMatrix_);

        if (conservative)
        {
            //-mass conservative terms
            hMEqn += (matrixPcModel_->Ch()*(hMatrix_.oldTime()-hMatrix_.prevIter())
                + ( Smatrix_ - Smatrix_.oldTime()))
                /hMatrix_.time().deltaT();
        }
    }

    
}

// ************************************************************************* //
