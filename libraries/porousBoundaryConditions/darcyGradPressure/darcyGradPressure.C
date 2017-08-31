/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "darcyGradPressure.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyGradPressure::darcyGradPressure
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_("Mf"),
    phiName_("phi"),
    phiGfName_("phiG"),
    phiPcName_("phiPc")
{}

Foam::darcyGradPressure::darcyGradPressure
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    phiGfName_(dict.lookupOrDefault<word>("phiG","phiG")),
    phiPcName_(dict.lookupOrDefault<word>("phiPc","phiPc"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}

Foam::darcyGradPressure::darcyGradPressure
(
    const darcyGradPressure& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_),
    phiGfName_(ptf.phiGfName_),
    phiPcName_(ptf.phiPcName_)
{}

Foam::darcyGradPressure::darcyGradPressure
(
    const darcyGradPressure& ptf
)
    :
    fixedGradientFvPatchScalarField(ptf),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_),
    phiGfName_(ptf.phiGfName_),
    phiPcName_(ptf.phiPcName_)
{}

Foam::darcyGradPressure::darcyGradPressure
(
    const darcyGradPressure& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_),
    phiGfName_(ptf.phiGfName_),
    phiPcName_(ptf.phiPcName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyGradPressure::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& Mf=
        patch().lookupPatchField<surfaceScalarField, scalar>(MfName_);

    const fvsPatchField<scalar>& phi=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvsPatchField<scalar>& phiGf=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiGfName_);

    const fvsPatchField<scalar>& phiPc=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiPcName_);

    //Extract the dictionary from database
    scalar  activateCapillarity(db().lookupObject<dictionary>("transportProperties").lookupOrDefault<scalar>("activateCapillarity",0.));

    gradient() = - (phi-phiGf-phiPc*activateCapillarity)/Mf/(patch().magSf());

    fixedGradientFvPatchScalarField::updateCoeffs();
}

void Foam::darcyGradPressure::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "phiG", "phiG", phiGfName_);
    writeEntryIfDifferent<word>(os, "phiPc", "phiPc", phiPcName_);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyGradPressure
);
}

// ************************************************************************* //
