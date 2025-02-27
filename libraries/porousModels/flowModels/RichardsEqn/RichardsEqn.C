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
#include "fvm.H"
#include "fvc.H"
#include "linear.H"
#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowModels::RichardsEqn::RichardsEqn
(
    const fvMesh& mesh,
    const IOdictionary& transportProperties,
    twophasePorousMediumModel& pmModel,
    incompressiblePhase& fluidPhase,
    const bool steady
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
    steady_(steady),
    rho_(fluidPhase.rho()),
    mu_(fluidPhase.mu()),
    Ss_(transportProperties.lookupOrDefault<dimensionedScalar>("Ss",dimensionedScalar("Ss",dimless/dimLength,0.))),
    phi_("phi", fluidPhase.phi()),
    Kf_(fvc::interpolate(K_,"K")),
    krf_("krthetaf",fvc::interpolate(krModel_.krb(),"krtheta")),
    Lf_("Lf",rho_*Kf_*krf_/mu_),
    Mf_("Mf",mag(g_)*Lf_),
    phiG_("phiG",(Lf_ * g_) & mesh.Sf()),
    phiPc_("phiPc",0*phiG_),
    patchDEM_(transportProperties.lookupOrDefault<word>("patchDEM","none")),
    patchDEMID_(mesh.boundaryMesh().findPatchID(patchDEM_)),
    topCellID_(0),
    seepageIDList_(0),
    distanceToDEM_(0),
    seepageValueList_(0)
{
    //- initialization
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

    //- Check Seepage condition
    if (patchDEM_ == "none")
    {
        Info << nl << "no DEM patch (no seepage condition)" << endl;
    }
    else
    {
        if (patchDEMID_ == -1)
        {
            FatalErrorIn("RichardsEqn.C") << "patch for seepage : " << patchDEM_ << " not found" << abort(FatalError);
        }
        else
        {
            Info << nl << "DEM patch used for seepage = " << patchDEM_ << " (id=" << patchDEMID_ << ")" << endl;
            topCellID_.resize(mesh.boundaryMesh()[patchDEMID_].size());
            topCellID_ = mesh.boundaryMesh()[patchDEMID_].faceCells();
            distanceToDEM_.resize(mesh.boundaryMesh()[patchDEMID_].size());
            distanceToDEM_ = mag(mesh.boundary()[patchDEMID_].delta()().component(2));
        }
    }
    updateProperties(true);
}
// * * * * * * * * * * * * * * * * * Members * * * * * * * * * * * * * * * * //

void Foam::flowModels::RichardsEqn::updateProperties(bool derivative)
{
    theta_ = pcModel_.correctAndSb(h_);
    krModel_.correctkrb(theta_, derivative);
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

void Foam::flowModels::RichardsEqn::updateSeepage()
{
    if (patchDEMID_ > -1)
    {
        seepageIDList_.clear();
        seepageValueList_.clear();
        volScalarField cellFlux(fvc::div(phi_));
        forAll(topCellID_,celli)
        {
            label currentCell = topCellID_[celli];

            if(h_[currentCell] >= distanceToDEM_[celli])
            {
                if (cellFlux[currentCell] < 0)
                {
                    seepageIDList_.append(currentCell);
                    seepageValueList_.append(distanceToDEM_[celli]);
                }
            }
        }

        // Display number of seepage cells
        label nSeepageCells = seepageIDList_.size();
        reduce(nSeepageCells, sumOp<label>());
        Info << "Number of seepage cells = " << nSeepageCells << endl;
    }
}

void Foam::flowModels::RichardsEqn::noConvergence
(
    multiDtManager& MDTM,
    Time& runTime,
    label algoID
)
{
    //- Set h equal to old-time value
    h_ = h_.oldTime();
    pmModel_.rewindTime();

    //- Rewind time
    runTime.setTime(runTime.timeOutputValue()-runTime.deltaTValue(),runTime.timeIndex());

    //- Reset iterator to indicate non-convergence and update timestep
    MDTM.dtManagerI(algoID).reset(MDTM.dtManagerI(algoID).maxIter()+1);
    MDTM.updateDt();

    //- Update time output value
    runTime.setTime(runTime.timeOutputValue()+runTime.deltaTValue(),runTime.timeIndex());

    //- Update properties
    updateProperties(true);
}

Foam::fvScalarMatrix Foam::flowModels::RichardsEqn::buildEqn()
{
    pmModel_.correct(h_, steady_, massConservative_);

    h_.storePrevIter();

    fvScalarMatrix hEqn
        (
            //- transport terms
            - fvm::laplacian(Mf_,h_)
            + fvc::div(phiG_)
            ==
            - pmModel_.exchangeTerm()
            - sourceTerm_
        );

    if (!steady_)
    {
        //- accumulation terms
        hEqn += (Ss_*pcModel_.Se() + pcModel_.Ch()) * fvm::ddt(h_);

        if (massConservative_)
        {
            //-mass conservative terms
            hEqn += (pcModel_.Ch()*(h_.oldTime()-h_.prevIter())
            + (theta_ - theta_.oldTime())) / h_.mesh().time().deltaT();
        }
    }

    if (seepageIDList_.size() > 0)
        hEqn.setValues(seepageIDList_,seepageValueList_);

    return hEqn;
}

Foam::scalar Foam::flowModels::RichardsEqn::initResidual
(
    const fvScalarMatrix& hEqn
) {
    scalarField Ax(mesh_.nCells());
    scalarField Axbar(mesh_.nCells());
    scalar xbar(gAverage(h_));

    const scalar *const __restrict__ diagPtr = hEqn.diag().begin();
    const label *const __restrict__ uPtr = hEqn.lduAddr().upperAddr().begin();
    const label *const __restrict__ lPtr = hEqn.lduAddr().lowerAddr().begin();
    const scalar *const __restrict__ upperPtr = hEqn.upper().begin();
    const label nCells = hEqn.diag().size();

    for (label cell = 0; cell < nCells; cell++) {
        Ax[cell] = diagPtr[cell] * h_[cell];
        Axbar[cell] = diagPtr[cell] * xbar;
    }
    const label nFaces = hEqn.upper().size();
    for (label face = 0; face < nFaces; face++) {
        Ax[uPtr[face]] += upperPtr[face] * h_[lPtr[face]];
        Ax[lPtr[face]] += upperPtr[face] * h_[uPtr[face]];
        Axbar[uPtr[face]] += upperPtr[face] * xbar;
        Axbar[lPtr[face]] += upperPtr[face] * xbar;
    }

    forAll(hEqn.internalCoeffs(), patchi)
    {
        forAll(hEqn.internalCoeffs()[patchi], facei)
        {
            label id_cell = h_.mesh().boundary()[patchi].faceCells()[facei];
            Ax[id_cell] += hEqn.internalCoeffs()[patchi][facei] * h_[id_cell];
            Ax[id_cell] -= hEqn.boundaryCoeffs()[patchi][facei];
            Axbar[id_cell] += hEqn.internalCoeffs()[patchi][facei] * xbar;
            Axbar[id_cell] -= hEqn.boundaryCoeffs()[patchi][facei];
        }
    }

    scalar normFactor = gSumMag(Ax - Axbar) + gSumMag(hEqn.source() - Axbar);
    return gSumMag(hEqn.source() - Ax) / normFactor;
}


const Foam::Tuple2<Foam::scalar, Foam::scalar> Foam::flowModels::RichardsEqn::solvePicard
(
    const scalar tolerance
)
{
    Tuple2<scalar, scalar> res(0, 0);
    fvScalarMatrix hEqnPicard = buildEqn();
    res.first() = hEqnPicard.solve().initialResidual();
    if (res.first() > tolerance) h_.relax();
    scalarField deltah(h_-h_.prevIter());
    forAll(seepageIDList_,celli) deltah[seepageIDList_[celli]] = 0;
    res.second() = gMax(mag(deltah)());
    return res;
}

const Foam::Tuple2<Foam::scalar, Foam::scalar> Foam::flowModels::RichardsEqn::solveNewton
(
    const scalar tolerance
)
{
    //- Compute initial residual
    Tuple2<scalar, scalar> res(0, 0);
    fvScalarMatrix hEqnPicard = buildEqn();
    res.first() = initResidual(hEqnPicard);
    tmp<DimensionedField<scalar, volMesh>> ResiduN = DimensionedField<scalar, volMesh>::New(
        "ResiduN",
        mesh_,
        dimless/dimTime,
        (hEqnPicard.residual()/mesh_.V())()
    );
    Info << "Richards' equation initial residual = " << res.first() << endl;

    //- Construct and solve Newton system
    volScalarField dLdS( pcModel_.Ch() * rho_ * K_ * krModel_.dkrbdS() / mu_ );
    volScalarField dMdS( mag(g_) * dLdS );

    fvScalarMatrix deltahEqn_hGrad(dMdS  * fvm::div(fvc::snGrad(h_) * mesh_.magSf(), deltah_));
    fvScalarMatrix deltahEqn_grav( dLdS * fvm::div(-g_ & mesh_.Sf(), deltah_));

    deltahEqn_hGrad.diag() *= -1;
    deltahEqn_grav.diag() *= -1;

    deltah_.storePrevIter();

    fvScalarMatrix deltahEqn
        (
            - fvm::laplacian(Mf_, deltah_)
            + deltahEqn_hGrad
            + deltahEqn_grav
            ==
            ResiduN
        );

    if (!steady_)
    {
        deltahEqn += fvm::Sp(
                (
                    (Ss_ * pcModel_.Ch()) * (h_ - h_.oldTime())
                    + Ss_ * pcModel_.Se()
                    + pcModel_.Ch()
                ) / mesh_.time().deltaT()
                , deltah_);
    }

    scalarField forInversion = deltahEqn.upper();
    deltahEqn.upper() = deltahEqn.lower();
    deltahEqn.lower() = forInversion;

    forAll(mesh_.boundary(),patchi)
    {
        if (deltah_.boundaryField().types()[patchi] == "darcyGradPressure")
        {
            deltahEqn.boundaryCoeffs()[patchi] = 0;
        }
    }
    if (seepageIDList_.size() > 0) {
        deltahEqn.setValues(seepageIDList_,0);
    }

    deltahEqn.solve();
    if (res.first() > tolerance) deltah_.relax();
    h_ = h_.prevIter() + deltah_;
    h_.correctBoundaryConditions();

    res.second() = gMax(mag(deltah_)());
    return res;
}


void Foam::flowModels::RichardsEqn::info()
{
    scalarField dtheta_tmp = mag(theta_.internalField()-theta_.oldTime().internalField());
    scalar dtheta = gMax(dtheta_tmp);

    //- water mass balance terminal display
    Info << nl << "Saturation theta: min(theta) = " << gMin(theta_.internalField())
                          << " max(theta) = " << gMax(theta_.internalField()) << " dthetamax = " << dtheta << endl;
    Info << "Head pressure h: min(h) = " << gMin(h_.internalField())
                         << " max(h) = " << gMax(h_.internalField()) << endl;
    Info << "Water mass balance (m3/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm_).value() << " ; ";
    forAll(phi_.boundaryField(),patchi)
    {
        if (mesh_.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi_.boundaryField()[patchi].patch().name() << " = " <<  gSum(phi_.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;
}

void Foam::flowModels::RichardsEqn::checkSteadyConfig
(
    const scalar tol1,
    const scalar tol2
)
{
    if (steady_ && tol1 < 1.0 && tol2 < 1.0)
    {
        FatalErrorIn("RichardsEqn.C") << "Only one tolerance (Newton or Picard) should be defined in fvSolution "
            << " for steady simulations" << abort(FatalError);
    }
};

// ************************************************************************* //
