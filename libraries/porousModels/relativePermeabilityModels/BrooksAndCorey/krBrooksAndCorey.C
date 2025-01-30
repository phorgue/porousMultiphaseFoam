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

#include "krBrooksAndCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krBrooksAndCorey, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krBrooksAndCorey,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krBrooksAndCorey::krBrooksAndCorey
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname,
    const word porousRegion
)
    :
    relativePermeabilityModel(mesh, transportProperties.subDict(typeName + "Coeffs"),
                              Sname, porousRegion),
    n_
    (
        IOobject
        (
            "n"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.getOrDefault<scalar>("n"+porousRegion,0))
    ),
    kramax_
    (
        IOobject
        (
            "kr"+Sname+"max"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.getOrDefault<scalar>("kr"+Sname+"max"+porousRegion,1.0))
    ),
    krbmax_
    (
        IOobject
        (
            "kr"+Sname+"max"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.getOrDefault<scalar>("kr"+Sname+"max"+porousRegion,1.0))
    )
{
    if (gMin(n_) <= 0)
    {
        FatalErrorIn
            (
                "in krBrooksAndCorey.C"
            )
            << "Relative permeability coefficient n equal or less than 0"
                << nl << "n" + porousRegion << " is required"
                << exit(FatalError);
    }
    Info << "Brooks and Corey parameters for relative permeability model" << nl << "{" << endl;
     Info <<  "    " << Sname << porousRegion << "min" << " ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    " << Sname << porousRegion << "max" << " ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
   Info << "    n" << porousRegion << " ";
    if (n_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(n_).value() << endl;}
    Info << "    kramax" << porousRegion << " ";
    if (kramax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(kramax_).value() << endl;}
    Info << "    krbmax" << porousRegion << " ";
    if (krbmax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(krbmax_).value() << endl;}
    Info << "} \n" << endl;
}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * //

void Foam::relativePermeabilityModels::krBrooksAndCorey::correct(const volScalarField& Sb, bool derivative)
{
    Se_= (Sb-Smin_)/(Smax_-Smin_);
    kra_ = kramax_* pow((scalar(1)-Se_),n_);
    krb_ = krbmax_ * pow(Se_,n_);
    if (derivative)
    {
        dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
        dkrbdS_ = krbmax_*n_*pow(Se_,n_-1)/(Smax_- Smin_);
    }
}
void Foam::relativePermeabilityModels::krBrooksAndCorey::correctkra(const volScalarField& Sb, bool derivative)
{
    Se_= (Sb-Smin_)/(Smax_-Smin_);
    kra_ = kramax_* pow((scalar(1)-Se_),n_);
    if (derivative)
    {
        dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);

    }
}
void Foam::relativePermeabilityModels::krBrooksAndCorey::correctkrb(const volScalarField& Sb, bool derivative)
{
    Se_= (Sb-Smin_)/(Smax_-Smin_);
    krb_ = krbmax_ * pow(Se_,n_);
    if (derivative)
    {
        dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
    }
}
void Foam::relativePermeabilityModels::krBrooksAndCorey::correctkrb(const volScalarField& Sb, const label& celli)
{
    scalar Se = (Sb[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
    krb_[celli] = krbmax_[celli] * pow(Se,n_[celli]);
}
Foam::tmp<Foam::volScalarField> Foam::relativePermeabilityModels::krBrooksAndCorey::kr(const volScalarField& S)
{
    volScalarField Se((S-Smin_)/(Smax_-Smin_));
    return krbmax_ * pow(Se,n_);
}

// ************************************************************************* //
