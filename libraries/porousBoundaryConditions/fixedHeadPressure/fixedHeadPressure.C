/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "fixedHeadPressure.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedHeadPressure::
fixedHeadPressure
(
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(h, iF),
    potential_(0.)
{}


Foam::fixedHeadPressure::
fixedHeadPressure
(
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(h, iF, dict, false),
    potential_(dict.lookupOrDefault<scalar>("potential",0.))
{}


Foam::fixedHeadPressure::
fixedHeadPressure
(
    const fixedHeadPressure& ptf,
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, h, iF, mapper),
    potential_(ptf.potential_)
{}


Foam::fixedHeadPressure::
fixedHeadPressure
(
    const fixedHeadPressure& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    potential_(ptf.potential_)
{}


Foam::fixedHeadPressure::
fixedHeadPressure
(
    const fixedHeadPressure& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    potential_(ptf.potential_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedHeadPressure::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const vectorField& fp = patch().patch().faceCentres();
    scalarField results(patch().patch().faceCentres().size());    
    forAll(fp,facei)
    {
        results[facei] = (potential_ - fp[facei].z());
    }
    
    operator== (results);
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fixedHeadPressure::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedHeadPressure
    );
}

// ************************************************************************* //
