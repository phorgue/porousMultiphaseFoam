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

#include "eventFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eventFlux::
eventFlux
(
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedValueFvPatchScalarField(h, iF),
    eventFluxValue_(0.),
    phiName_("phi"),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{}


Foam::eventFlux::
eventFlux
(
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedValueFvPatchScalarField(h, iF, dict, false),
    eventFluxValue_(dict.lookupOrDefault<scalar>("constantValue",0.)),
    phiName_(dict.lookupOrDefault<word>("phiName","phi")),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{
    word eventFileName = dict.lookupOrDefault<word>("eventFile","");
    //- Read if backward time scheme is used
    if (word(internalField().mesh().ddtScheme("source")) == "backward")
    {
        isBackwardScheme_ = true;
    }
    
    if (eventFileName != "")
    {
        //- reading patch event file, compute current value, store to old values
        eventFile_.read(eventFileName,true);
        eventFile_.update(this->db().time().startTime().value());
        eventFile_.storeOldValues();

        //- Reading patch event file and adding intermediate time step
        scalar eventTimeStep = this->db().time().controlDict().lookupOrDefault<scalar>("eventTimeStep",0);
        if (eventTimeStep > 0)
        {
            eventFile_.addIntermediateTimeSteps(eventTimeStep);
        }

        //- finding corresponding patch name
        forAll(eventFile_.patchNameList(),patchi)
        {
            if (patch().name() == eventFile_.patchNameList()[patchi]) patchEventID_ = patchi;
        }
        if (patchEventID_ == -1)
        {
            FatalErrorIn("eventFlux.C") << " patch '" << patch().name() << "' not found in event file : " << eventFile_.name() << abort(FatalError);
        }
    }
    else
    {
        Info << "eventFlux boundary condition without event file" << endl;
    }

    //- Read if backward time scheme is used
    if (word(internalField().mesh().ddtScheme("source")) == "backward")
    {
        isBackwardScheme_ = true;
    }

}


Foam::eventFlux::
eventFlux
(
    const eventFlux& ptf,
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, h, iF, mapper),
    eventFluxValue_(ptf.eventFluxValue_),
    phiName_(ptf.phiName_),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{}


Foam::eventFlux::
eventFlux
(
    const eventFlux& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    eventFluxValue_(ptf.eventFluxValue_),
    phiName_(ptf.phiName_),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{}


Foam::eventFlux::
eventFlux
(
    const eventFlux& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    eventFluxValue_(ptf.eventFluxValue_),
    phiName_(ptf.phiName_),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eventFlux::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const fvsPatchField<scalar>& phip_=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    scalar valueEvent = 0.0;

    if (patchEventID_ != -1)
    {
        if (isBackwardScheme_)
        {
            scalar deltaT =this->db().time().deltaT().value();
            scalar deltaT0 =this->db().time().deltaT0().value();
            scalar coefft0_00 = deltaT/(deltaT + deltaT0);
            scalar coefftn_0 = 1 + coefft0_00;
            valueEvent = coefftn_0*eventFile_.currentValue(patchEventID_) - coefft0_00*eventFile_.oldValue(patchEventID_);
        }
        else
        {
            valueEvent = eventFile_.currentValue(patchEventID_);
        }

        //- Updating event value
        eventFile_.update(this->db().time().value());
    }

    //- Computing fixed value
    scalarField results(patch().patch().faceCentres().size());    
    results = (eventFluxValue_+valueEvent)/sum(phip_);
    operator== (results);
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::eventFlux::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    eventFlux
);
}

// ************************************************************************* //
