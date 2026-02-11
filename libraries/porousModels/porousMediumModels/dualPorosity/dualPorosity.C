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
    twophasePorousMediumModel,
    dualPorosity,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumModels::dualPorosity::dualPorosity
(
    const word& Sname,
    const fvMesh& mesh,
    const IOdictionary& transportProperties,
    const autoPtr<incompressiblePhase>& phase,
    const word& porousRegion
)
    :
    twophasePorousMediumModel(Sname, mesh, transportProperties, phase, "Fracture"),
    dualPorosityCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    mesh_(mesh),
    phase_(phase),
    g
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),
    hMatrix_
    (
        IOobject
        (
            "hMatrix",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Smatrix_
    (
        IOobject
        (
            Sname_+"Matrix",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless,
        calculatedFvPatchScalarField::typeName
    ),
    KMatrix_
    (
        IOobject
        (
            "KMatrix",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dualPorosityCoeffs_.lookupOrDefault("KMatrix", dimensionedScalar(dimArea, 0))
    ),
    KMatrixf_(fvc::interpolate(KMatrix_, "K")),
    KExchange_
    (
        IOobject
        (
            "KExchange",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dualPorosityCoeffs_.lookupOrDefault("KExchange", dimensionedScalar(dimArea, 0))
    ),
    a_(dualPorosityCoeffs_.get<dimensionedScalar>("a")),
    beta_(dualPorosityCoeffs_.get<dimensionedScalar>("beta")),
    gammaW_(dualPorosityCoeffs_.get<dimensionedScalar>("gammaW")),
    geomFactor_(beta_/(a_*a_)*gammaW_),
    UMatrix_
    (
        IOobject
        (
            "U"+Sname_+"Matrix",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    LMatrixf_(phase_->rho()*KMatrixf_/phase_->mu()),
    MMatrixf_(mag(g)*LMatrixf_),
    phiGMatrixf_((LMatrixf_ * g) & mesh_.Sf()),
    phiMatrix_
    (
        IOobject
        (
            "phiMatrix",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(UMatrix_) & mesh.Sf()
    )
{
    Info << "Dual porosity parameters for water flow " << nl << "{" << endl;
    Info << "    a " << a_.value() << endl;
    Info << "    beta " << beta_.value() << endl;
    Info << "    gammaW " << gammaW_.value() << endl;
    check_K(KMatrix_, dualPorosityCoeffs_);
    check_K(KExchange_, dualPorosityCoeffs_);
    Info << "} \n" << endl;
    sourceTerm_.writeOpt(IOobject::AUTO_WRITE);
    matrixPcModel_ = capillarityModel::New(mesh, transportProperties, Sname_, "Matrix");
    matrixKrModel_ = relativePermeabilityModel::New(mesh, transportProperties, Sname_, "Matrix");
    updateMatrixProperties();
    Smatrix_.write();
}

// * * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

void Foam::porousMediumModels::dualPorosity::updateMatrixProperties()
{
    //- update matrix properties
    Smatrix_ = matrixPcModel_->correctAndSb(hMatrix_);
    matrixKrModel_->correctkrb(Smatrix_, false);
    surfaceScalarField krthetaMatrixf(fvc::interpolate(matrixKrModel_->krb(), "krthetaMatrix"));
    LMatrixf_ = phase_->rho()*KMatrixf_*krthetaMatrixf/phase_->mu();
    MMatrixf_ = mag(g)*LMatrixf_;
    phiGMatrixf_ = (LMatrixf_ * g) & mesh_.Sf();
    phiMatrix_ = phiGMatrixf_-(MMatrixf_*fvc::snGrad(hMatrix_))*mesh_.magSf();
    UMatrix_ = fvc::reconstruct(phiMatrix_);
    UMatrix_.correctBoundaryConditions();
    forAll(UMatrix_.boundaryField(),patchi)
    {
        if (isA< fixedValueFvPatchField<Vector<scalar>> >(UMatrix_.boundaryField()[patchi]))
        {
            phiMatrix_.boundaryFieldRef()[patchi] = UMatrix_.boundaryField()[patchi] & mesh_.Sf().boundaryField()[patchi];
        }
    }
}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

void Foam::porousMediumModels::dualPorosity::rewindTime()
{
    hMatrix_ = hMatrix_.oldTime();
}

void Foam::porousMediumModels::dualPorosity::correct()
{
    FatalErrorIn("dualPorosity.C") << " dualPorosity cannot be used with impesFoam/anisoImpesFoam " << abort(FatalError);
}

void Foam::porousMediumModels::dualPorosity::correct(volScalarField& hFracture, const bool steady, const bool conservative)
{
    hMatrix_.storePrevIter();

    //- Update transfert coefficient
    volScalarField SFracture(matrixPcModel_->S(hFracture));
    volScalarField krExchange((matrixKrModel_->kr(SFracture)+matrixKrModel_->krb())/2.0);
    volScalarField alphaW(geomFactor_*phase_->rho()*mag(g)*KExchange_*krExchange/phase_->mu());

    //- solve matrix equation
    fvScalarMatrix hMEqn
        (
            //- transport terms
            - fvm::laplacian(MMatrixf_,hMatrix_)
            + fvc::div(phiGMatrixf_)
            ==
            alphaW * hFracture
            - fvm::Sp(alphaW, hMatrix_)
        );

    if (!steady)
    {
        //- accumulation terms
        hMEqn += matrixPcModel_->Ch() * fvm::ddt(hMatrix_);

        if (conservative)
        {
            //-mass conservative terms
            hMEqn += (matrixPcModel_->Ch()*(hMatrix_.oldTime()-hMatrix_.prevIter())
                + ( Smatrix_ - Smatrix_.oldTime())) / hMatrix_.time().deltaT();
        }
    }

    if (seepageIDList_->size() > 0) hMEqn.setValues(*seepageIDList_,*seepageValueList_);

    hMEqn.solve();

    //- compute source term using update hMatrix field
    exchangeTerm_ = alphaW* (hFracture - hMatrix_);

    //- update properties using new solution
    updateMatrixProperties();
}

// ************************************************************************* //
