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

#include "DupuitDarcyEqn.H"
#include "fvm.H"
#include "fvc.H"
#include "linear.H"
#include "fixedValueFvPatchField.H"
#include "DEMfile.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowModels::DupuitDarcyEqn::DupuitDarcyEqn
(
    const fvMesh& mesh,
    const IOdictionary& transportProperties,
    porousMediumModel& pmModel,
    incompressiblePhase& fluidPhase,
    const bool steady
)
    :
    g_("g",dimLength/(dimTime*dimTime),9.81),
    potential_
    (
        IOobject
        (
            "potential",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    z0_
    (
        IOobject
            (
                "z0",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
        mesh
    ),
    DEM_
    (
        IOobject
        (
            "potentialDEM",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        z0_
    ),
    hwater_
    (
        IOobject
            (
                "hwater",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
        potential_ - z0_
    ),
    infiltration_
    (
        IOobject
        (
            "infiltration",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        transportProperties.lookupOrDefault("infiltration", dimensionedScalar("", dimLength/dimTime, 0.))
    ),
    seepageTerm_
    (
        IOobject
        (
            "seepageTerm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("seepage_value",dimLength/dimTime,0.)
    ),
    mesh_(mesh),
    pmModel_(pmModel),
    sourceTerm_(pmModel_.sourceTerm()),
    eps_(pmModel_.eps()),
    K_(pmModel_.K()),
    U_(fluidPhase.U()),
    steady_(steady),
    seepage_(transportProperties.lookupOrDefault("seepage", false)),
    DEMfileName_(transportProperties.lookupOrDefault<word>("fileDEM","")),
    zScale_(mesh.bounds().max().z()-mesh.bounds().min().z()),
    rho_(fluidPhase.rho()),
    mu_(fluidPhase.mu()),
    hwaterMin_(transportProperties.lookupOrDefault("hwaterMin",dimensionedScalar("hwaterMin",dimLength,0.01))),
    cumulativeWaterAdded_(0),
    flowInOutFixedPoints_(0),
    flowOutSeepage_(0),
    phi_("phi", fluidPhase.phi()),
    Kf_(fvc::interpolate(K_,"K")),
    Mf_("Mf",Kf_*g_*rho_/mu_),
    transmissivity_("transmissivity",Mf_*fvc::interpolate(hwater_)),
    phiG_("phiG", 0*phi_),
    phiPc_("phiPc",0*phi_),
    phihwater_(phi_ * fvc::interpolate(hwater_)),
    seepageIDList_(0),
    dryCellIDList_(0),
    fixedPotentialIDList_(0),
    seepageValueList_(0),
    fixedPotentialValueList_(0),
    cellFlux_(fvc::div(phi_ * fvc::interpolate(hwater_)))
{
    //- Display initial information
    Info << nl << "Elevation field: min(z0) = " << gMin(z0_.internalField())
         << " ; max(z0) = " << gMax(z0_.internalField()) << endl;

    //- Check hwater/potential values
    correctInitialPotential();
    Info << nl << "Water depth: min(hwater) = " << gMin(hwater_.internalField())
         << " ; max(hwater) = " << gMax(hwater_.internalField()) << endl;

    if (!hwater_.headerOk()) hwater_.write();

    //- Check porosity/permeability
    if (!steady) pmModel_.check_eps();
    pmModel_.check_K();

    //- Computing DEM if necessary
    if (DEM_.headerOk()) Info << nl << "Reading precomputed potentialDEM file in constant/" << endl;
    else
    {
        if (DEMfileName_.size() > 0)
        {
            Info << nl << "Reading DEM file to compute potentialDEM...";
            DEMfile potentialDEMfile(DEMfileName_);
            Info << "OK" << endl;

            //- Computing for potentialDEM for Seep term
            Info << "Interpolating value for potentialDEM...";
            forAll(mesh.C(),celli)
            {
                DEM_[celli] = potentialDEMfile.interpolate(mesh.C()[celli]);
            }
            Info << "OK" << endl;
            DEM_.write();
        }
        else Info << nl << "no potentialDEM field";
    }

    //- checking that potentialDEM in superior to z0)
    if (gMin((DEM_.internalField() - z0_.internalField())()) < 0)
    {
        Warning() << "potential DEM inferior to z0 in domain => set to z0+hwaterMin" << endl;
        forAll(DEM_.internalField(),celli)
        {
            DEM_.ref()[celli] = max(DEM_.internalField()[celli],z0_.internalField()[celli]+hwaterMin_.value());
        }
    }

    //- Checking seepage option
    if (seepage_)
    {
        Info << nl << "Seepage option is active" << endl;
        if (!DEM_.headerOk() && DEMfileName_.size() == 0)
        {
            FatalErrorIn("readFixedPoints.H") << nl << "no potentialDEM file neither DEM file while seepage is active " << abort(FatalError);
        }
    }
    else Info << nl << "Seepage is inactive" << endl;

    //- reading fixed potential list if present
    List<Tuple2<point,scalar> > fixedPotentialList(transportProperties.lookupOrDefault("fixedPotentialList",List<Tuple2<point,scalar> >()));
    bool useDEMtoFixPotential(transportProperties.lookupOrDefault<bool>("useDEMtoFixPotential",false));
    if (fixedPotentialList.size() > 0) {
        initFixedPotential(fixedPotentialList, useDEMtoFixPotential);
    }

    updateProperties();
}
// * * * * * * * * * * * * * * * * * Members * * * * * * * * * * * * * * * * //

void Foam::flowModels::DupuitDarcyEqn::initFixedPotential(const List<Tuple2<point,scalar>>& fpl, const bool useDEM)
{
    //- find the closest cell to fixed points
    label nFixedPoints = 0;
    forAll(fpl, pointi)
    {
        label cellID =  mesh_.findCell(fpl[pointi].first());
        if (cellID > -1)
        {
            nFixedPoints += 1;
            fixedPotentialIDList_.resize(nFixedPoints, cellID);
            fixedPotentialValueList_.resize(nFixedPoints, fpl[pointi].second());
        }
        else
        {
            Warning() << "fixed point " << fpl[pointi].first() << " is not in the domain";
        }
    }

    if (useDEM)
    {
        Info << nl << "potentialDEM used for fixedPotentialList" << endl;
        forAll(fpl, pointi)
        {
            fixedPotentialValueList_[pointi] = DEM_[fixedPotentialIDList_[pointi]];
        }
    }
    else
    {
        Info << nl << "user-defined values used for fixedPotentialList" << endl;
    }

    //- Display information about fixed values
    if (Pstream::nProcs() == 1)
    {
        Info << nl << "Fixed potential positions and values are " << nl << "{";
        forAll(fpl,pointi)
        {
            scalar distance_to_centre = Foam::sqrt(pow(mesh_.C()[fixedPotentialIDList_[pointi]].x()-fpl[pointi].first().x(),2)
                                                   +pow(mesh_.C()[fixedPotentialIDList_[pointi]].y()-fpl[pointi].first().y(),2));
            Info << nl << "  " << fpl[pointi].first() <<  " : value " << fixedPotentialValueList_[pointi]
                 << "  (cellID = " << fixedPotentialIDList_[pointi] << ", distance with cell-center = " << distance_to_centre << ")";
        }
        Info << nl << "}" << endl;
    }
}

void Foam::flowModels::DupuitDarcyEqn::correctInitialPotential()
{
    if (gMin(hwater_.internalField()) <= 0)
    {
        //- if cell initialization negative
        if (gMin(hwater_.internalField()) <= hwaterMin_.value())
        {
            Warning() << "Correcting hwater/potential cell values (when hwater < hwaterMin)" << endl;
            forAll(mesh_.C(),celli)
            {
                if ( hwater_[celli] < hwaterMin_.value())
                {
                    hwater_[celli] = hwaterMin_.value();
                    potential_[celli] = z0_[celli] + hwaterMin_.value();
                }
            }
        }

        //- if fixed boundary patch has negative hwater values
        forAll(mesh_.boundary(),patchi)
        {
            if (gMin(hwater_.boundaryField()[patchi]) <= hwaterMin_.value())
            {
                if (isA< fixedValueFvPatchField<scalar> >(hwater_.boundaryField()[patchi]))
                {
                    Warning() <<  "Correcting hwater/potential cell values (when hwater < hwaterMin) for patch "
                        << mesh_.boundary()[patchi].name() << endl;
                    forAll(hwater_.boundaryField()[patchi],facei)
                    {
                        if (hwater_.boundaryField()[patchi][facei] < hwaterMin_.value())
                        {
                             hwater_.boundaryFieldRef()[patchi][facei] = hwaterMin_.value();
                             potential_.boundaryFieldRef()[patchi][facei] =
                                 z0_.boundaryField()[patchi][facei] + hwaterMin_.value();
                        }
                    }
                }
            }
        }
    }
}

void Foam::flowModels::DupuitDarcyEqn::updateProperties()
{
    //- updating flow properties
    transmissivity_ = Mf_*fvc::interpolate(hwater_);
    phi_ = (-Mf_ * fvc::snGrad(potential_)) * mesh_.magSf();
    forAll(mesh_.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(U_.boundaryField()[patchi]))
        {
            phi_.boundaryFieldRef()[patchi] = U_.boundaryField()[patchi] & mesh_.Sf().boundaryField()[patchi];
        }
    }
    U_ = fvc::reconstruct(phi_);
    forAll(dryCellIDList_, celli)
    {
        U_[dryCellIDList_[celli]] = vector(0,0,0);
    }
    U_.correctBoundaryConditions();
    phihwater_ = phi_ * fvc::interpolate(hwater_);
    cellFlux_ = fvc::div(phihwater_) + infiltration_ + zScale_ * sourceTerm_;

    //- Compute outflow and seepage terms
    seepageTerm_.primitiveFieldRef() = 0;
    flowInOutFixedPoints_ = 0;
    if (fixedPotentialIDList_.size() > 0)
    {
        forAll(fixedPotentialIDList_, pointi)
        {
            label currentCell = fixedPotentialIDList_[pointi];
            scalar area = mesh_.V()[currentCell]/zScale_;
            if (cellFlux_.internalField()[currentCell] < 0) {
                seepageTerm_[currentCell] = - cellFlux_.internalField()[currentCell];
            }
            flowInOutFixedPoints_ -= cellFlux_.internalField()[currentCell]*area;
        }
    }
    flowOutSeepage_ = 0;
    if (seepage_)
    {
        seepageTerm_.primitiveFieldRef() = 0;
        forAll(seepageIDList_, pointi)
        {
            label currentCell = seepageIDList_[pointi];
            scalar area = mesh_.V()[currentCell]/zScale_;
            seepageTerm_[currentCell] = - cellFlux_.internalField()[currentCell];
            flowOutSeepage_ -= cellFlux_.internalField()[currentCell]*area;
        }
    }

}

void Foam::flowModels::DupuitDarcyEqn::updateSeepage()
{
    seepageIDList_.clear();
    seepageValueList_.clear();
    forAll(mesh_.C(),celli)
    {
        if(potential_[celli] >= DEM_[celli] && cellFlux_[celli] < 0)
        {
            seepageIDList_.append(celli);
            seepageValueList_.append(DEM_[celli]);
        }
    }

    // Display number of seepage cells
    label nSeepageCells = seepageIDList_.size();
    reduce(nSeepageCells, sumOp<label>());
    Info << "Number of seepage cells = " << nSeepageCells << endl;
}

void Foam::flowModels::DupuitDarcyEqn::updateDryCells()
{
    //- Checking for dry cells
    dryCellIDList_.clear();
    if (gMin(hwater_.internalField()) <= hwaterMin_.value())
    {
        scalar waterAdded = 0;
        forAll(hwater_, celli)
        {
            if (hwater_[celli] <= hwaterMin_.value())
            {
                dryCellIDList_.append(celli);
                waterAdded += (hwaterMin_.value()-hwater_[celli])*mesh_.V()[celli]/zScale_;
                hwater_[celli] = hwaterMin_.value();
            }
        }
        Info << "Number of dry cells = " << dryCellIDList_.size() << ", water added = " << waterAdded << " m3";
        if (!steady_)
        {
            cumulativeWaterAdded_ += waterAdded;
            Info << ", cumulative water added = " << cumulativeWaterAdded_ << " m3 ("
                 << cumulativeWaterAdded_ * zScale_ / gSum(mesh_.V()) << " m)" << endl;
        }
        Info << endl;
    }

}

Foam::scalar Foam::flowModels::DupuitDarcyEqn::solve ()
{
    potential_.storePrevIter();

    fvScalarMatrix potentialEqn
        (
            - fvm::laplacian(transmissivity_ ,potential_ ,"laplacian(transmissivity,potential)")
            ==
            - sourceTerm_ * zScale_
            - infiltration_
        );

    if (!steady_) potentialEqn += eps_ * fvm::ddt(potential_);

    //- Fixed potential values
    if (fixedPotentialIDList_.size() > 0)
    {
        potentialEqn.setValues(fixedPotentialIDList_,fixedPotentialValueList_);
    }

    //- Seepage
    if (seepage_)
    {
        updateSeepage();
        if (seepageIDList_.size() > 0) potentialEqn.setValues(seepageIDList_, seepageValueList_);
    }

    //- Solve equation
    scalar maxResidual = potentialEqn.solve().initialResidual();

    if (steady_) potential_.relax();

    //- Update properties + dry cells
    hwater_ = potential_ - z0_;
    updateDryCells();
    updateProperties();

    return maxResidual;
}


void Foam::flowModels::DupuitDarcyEqn::info()
{
    Info << "Potential min: " << gMin(potential_.internalField()) << ", max = " << gMax(potential_.internalField())
      << ", delta(potential) = "
      << gMax(mag(potential_.internalField()-potential_.oldTime().internalField())()) << endl;
}

// ************************************************************************* //
