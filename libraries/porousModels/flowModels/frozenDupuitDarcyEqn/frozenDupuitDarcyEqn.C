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

#include "frozenDupuitDarcyEqn.H"
#include "fvc.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowModels::frozenDupuitDarcyEqn::frozenDupuitDarcyEqn
(
    const fvMesh& mesh,
    const IOdictionary& transportProperties,
    porousMediumModel& pmModel,
    fluidPhase& fluid
)
    :
    hwater_
    (
        IOobject
        (
            "hwater",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("hwater_default",dimLength,0.)
    ),
    seepageTerm_
    (
        IOobject
        (
            "seepageTerm",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("seepage_value",dimLength/dimTime,0.)
    ),
    zScale_(mesh.bounds().max().z()-mesh.bounds().min().z()),
    phihwater_(fluid.phi() * fvc::interpolate(hwater_))
{

    //- Reading or re-constructing hwater field
    if (hwater_.headerOk())
    {
        Info << "min(hwater) = " << gMin(hwater_.internalField()) << " ; max(hwater) = " << gMax(hwater_.internalField()) << endl;
    }
    else
    {
        Info << "file hwater not found " << endl
             << nl << "Reading potential field..." << endl;
        volScalarField potential
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
        );

        Info << nl << "Reading field z0" << endl;
        volScalarField z0
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
        );
        Info << "min(z0) = " << gMin(z0.internalField()) << " ; max(z0) = " << gMax(z0.internalField()) << endl;

        hwater_ = potential - z0;

        Info << nl << "Computed hwater : min(hwater) = " << gMin(hwater_.internalField())
             << " ; max(hwater) = " << gMax(hwater_.internalField()) << endl;
        hwater_.write();
    }

    if (gMin(hwater_.internalField()) < 0)
    {
        FatalErrorIn("createFields") << " hwater has negative values" << abort(FatalError);
    }

    //- Reading or reconstructing phi
    volVectorField& U = fluid.U();
    surfaceScalarField& phi = fluid.phi();
    if (!exists(phi.objectPath()))
    {
        bool phiReconstruction = transportProperties.getOrDefault<bool>("phiReconstruction", true);
        Info << nl << "phi is not present: ";
        if (exists(mesh.time().timeName()/"potential") && exists(mesh.time().constant()/"K") &&
            transportProperties.found("rho") && transportProperties.found("mu") && phiReconstruction)
        {
            Info << "reconstructed using potential/permeability and fluid properties" << endl;
            volScalarField potential
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
            );
            volScalarField K
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
            );
            surfaceScalarField Kf = fvc::interpolate(K);
            dimensionedScalar g("g",dimLength/(dimTime*dimTime),9.81);
            const dimensionedScalar rho("rho",transportProperties);
            const dimensionedScalar mu("mu",transportProperties);
            phi = (-Kf*(g*rho/mu) * fvc::snGrad(potential)) * mesh.magSf();
            phi.write();
        }
        else
        {
            phi = Foam::linearInterpolate(U) & mesh.Sf();
            Info << "computed using linear interpolation from U (*WARNING* should not be conservative)" << endl;
            phi.write();
        }

        //- Recompute phihwater
        phihwater_ = phi * fvc::interpolate(hwater_);
    }


    if (seepageTerm_.headerOk())
    {
        Info << "seepageTerm field OK" << endl;
        //- ensuring that there is no seepage inflow
        forAll(seepageTerm_, celli) seepageTerm_[celli] = max(0, seepageTerm_[celli]);
    }
    else
    {
        Info << " no seepageTerm file" << endl;
    }

}

// * * * * * * * * * * * * * * * * * Members * * * * * * * * * * * * * * * * //

void Foam::flowModels::frozenDupuitDarcyEqn::info()
{
}

// ************************************************************************* //
