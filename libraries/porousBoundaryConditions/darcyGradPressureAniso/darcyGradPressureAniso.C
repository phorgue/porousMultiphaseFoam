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

#include "darcyGradPressureAniso.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "gravityMeshObject.H"
#include "linear.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_("Mf"),
    MbfName_("Mbf"),
    phiName_("phi"),
    LfName_("Lf"),
    gradpcName_("gradpc")
{}

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_(dict.getOrDefault<word>("Mf", "Mf")),
    MbfName_(dict.getOrDefault<word>("Mbf", "Mbf")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    LfName_(dict.getOrDefault<word>("Lf","Lf")),
    gradpcName_(dict.getOrDefault<word>("gradpc","gradpc"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const darcyGradPressureAniso& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    MfName_(ptf.MfName_),
    MbfName_(ptf.MbfName_),
    phiName_(ptf.phiName_),
    LfName_(ptf.LfName_),
    gradpcName_(ptf.gradpcName_)
{}

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const darcyGradPressureAniso& ptf
)
    :
    fixedGradientFvPatchScalarField(ptf),
    MfName_(ptf.MfName_),
    MbfName_(ptf.MbfName_),
    phiName_(ptf.phiName_),
    LfName_(ptf.LfName_),
    gradpcName_(ptf.gradpcName_)
{}

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const darcyGradPressureAniso& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    MbfName_(ptf.MbfName_),
    phiName_(ptf.phiName_),
    LfName_(ptf.LfName_),
    gradpcName_(ptf.gradpcName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyGradPressureAniso::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<tensor>& Mf=
        patch().lookupPatchField<surfaceTensorField, tensor>(MfName_);

    const fvsPatchField<tensor>& Mbf=
        patch().lookupPatchField<surfaceTensorField, tensor>(MbfName_);

    const fvsPatchField<scalar>& phi=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvsPatchField<tensor>& Lf=
        patch().lookupPatchField<surfaceTensorField, tensor>(LfName_);

    const fvsPatchField<vector>& gradpc=
        patch().lookupPatchField<surfaceVectorField, vector>(gradpcName_);

    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(db().time());

    //Activate or not capillarity term
    scalar  activateCapillarity(db().lookupObject<dictionary>("transportProperties").getOrDefault<scalar>("activateCapillarity",0.));
	
    gradient() = - (phi/patch().magSf()) * ( inv(Mf) & patch().nf() & patch().nf());
    gradient() = gradient() + (inv(Mf) & Lf & g.value() & patch().nf() );
    gradient() = gradient() + activateCapillarity*(inv(Mf) & Mbf & gradpc & patch().nf() );

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::darcyGradPressureAniso::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("Mf", "Mf", MfName_);
    os.writeEntryIfDifferent<word>("Mbf", "Mbf", MbfName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("Lf", "Lf", LfName_);  
    os.writeEntryIfDifferent<word>("gradpc", "gradpc", gradpcName_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyGradPressureAniso
);
}


// ************************************************************************* //
