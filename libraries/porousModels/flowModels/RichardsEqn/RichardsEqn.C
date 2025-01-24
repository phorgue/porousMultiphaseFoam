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

#include "RichardsEqn.H"
#include "fvc.H"
#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowModels::RichardsEqn::RichardsEqn
(
    const fvMesh& mesh,
    const IOdictionary& transportProperties,
    twophasePorousMediumModel& pmModel,
    incompressiblePhase& fluidPhase
)
    :
    g_(meshObjects::gravity::New(mesh.time())),
    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    deltah_
    (
        IOobject
        (
            "deltah",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        h_
    ),
    theta_
    (
        IOobject
        (
            "theta",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless,
        calculatedFvPatchScalarField::typeName
    ),
    mesh_(mesh),
    pmModel_(pmModel),
    krModel_(pmModel.krModel().ref()),
    pcModel_(pmModel.pcModel().ref()),
    sourceTerm_(pmModel_.sourceTerm()),
    K_(pmModel_.K()),
    U_(fluidPhase.U()),
    massConservative_(transportProperties.lookupOrDefault<bool>("massConservative",true)),
    rho_(fluidPhase.rho()),
    mu_(fluidPhase.mu()),
    Ss_(transportProperties.lookupOrDefault<dimensionedScalar>("Ss",dimensionedScalar("Ss",dimless/dimLength,0.))),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_) & mesh.Sf()
    ),
    Kf_(fvc::interpolate(K_,"K")),
    krf_("krthetaf",fvc::interpolate(krModel_.krb(),"krtheta")),
    Lf_("Lf",rho_*Kf_*krf_/mu_),
    Mf_("Mf",mag(g_)*Lf_),
    phiG_("phiG",(Lf_ * g_) & mesh.Sf()),
    phiPc_("phiPc",0*phiG_)
{
    //- initialization
    fluidPhase.phi().writeOpt()=IOobject::NO_WRITE;
    deltah_ == dimensionedScalar("",dimLength,0);

    //- Checking permeability field
    pmModel_.check_K();

    //- Checking gravity
    if (mag(g_).value() == 0)
    {
        FatalErrorIn("RichardsEqn.C")
            << " Magnitude of gravity mag(g) equal to zero " << abort(FatalError);
    }

    Info << nl << "Computing saturation field theta" << endl;
    theta_ = pcModel_.correctAndSb(h_);
    theta_.write();
}
// * * * * * * * * * * * * * * * * * Members * * * * * * * * * * * * * * * * //

void Foam::flowModels::RichardsEqn::updateProperties()
{
    theta_ = pcModel_.correctAndSb(h_);
    krModel_.correctkrb(theta_);
    krf_ = fvc::interpolate( krModel_.krb(),"krtheta");
    Lf_ = rho_ * Kf_ * krf_ / mu_;
    Mf_ = mag(g_) * Lf_;
    phiG_ = (Lf_ * g_) & mesh_.Sf();
    phi_ = phiG_ - (Mf_ * fvc::snGrad(h_)) * mesh_.magSf();
    U_ = fvc::reconstruct(phi_);
    U_.correctBoundaryConditions();
    forAll(mesh_.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(U_.boundaryField()[patchi]))
        {
            phi_.boundaryFieldRef()[patchi] = U_.boundaryField()[patchi] & mesh_.Sf().boundaryField()[patchi];
        }
    }
    h_.correctBoundaryConditions();
}

// ************************************************************************* //
