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

#include "fixedFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFlux::
fixedFlux
(
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedValueFvPatchScalarField(h, iF),
    fixedFluxValue_(0.),
    phiName_("phi"),
    iterEvent_(-1),
    valueEvent_(0.0)
{}


Foam::fixedFlux::
fixedFlux
(
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedValueFvPatchScalarField(h, iF, dict, false),
    fixedFluxValue_(dict.lookupOrDefault<scalar>("fixedFluxValue",0.)),
    phiName_(dict.lookupOrDefault<word>("phiName","phi")),
    iterEvent_(-1),
    valueEvent_(0.0)
{
    if (dict.found("event"))
    {
        List<Tuple2<scalar, scalar> > eventRead= dict.lookup("event");
        eventData_ = eventRead;
        iterEvent_ = 0;
    }
}


Foam::fixedFlux::
fixedFlux
(
    const fixedFlux& ptf,
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, h, iF, mapper),
    fixedFluxValue_(ptf.fixedFluxValue_),
    phiName_(ptf.phiName_),
    iterEvent_(-1),
    valueEvent_(0.0)
{}


Foam::fixedFlux::
fixedFlux
(
    const fixedFlux& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    fixedFluxValue_(ptf.fixedFluxValue_),
    phiName_(ptf.phiName_),
    iterEvent_(-1),
    valueEvent_(0.0)
{}


Foam::fixedFlux::
fixedFlux
(
    const fixedFlux& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    fixedFluxValue_(ptf.fixedFluxValue_),
    phiName_(ptf.phiName_),
    iterEvent_(-1),
    valueEvent_(0.0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFlux::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const fvsPatchField<scalar>& phip_=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    //- Updating event value
    if ((iterEvent_ >= 0 ) && (iterEvent_<eventData_.size()))
    {
        if (eventData_[iterEvent_+1].first() < this->db().time().value())
        {
            iterEvent_++;
            valueEvent_ = eventData_[iterEvent_].second();
        }
    }

    //- Computing fixed value
    scalarField results(patch().patch().faceCentres().size());    
    results = (-fixedFluxValue_-valueEvent_)/sum(phip_);
    operator== (results);
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fixedFlux::write(Ostream& os) const
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
        fixedFlux
    );
}

// ************************************************************************* //
