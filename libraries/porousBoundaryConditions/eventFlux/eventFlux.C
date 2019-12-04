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


Foam::List<Foam::patchEventFile*>* Foam::eventFlux::eventFileRegistry_ = nullptr;
Foam::word Foam::eventFlux::dtFieldNameOverride_ = "";

void Foam::eventFlux::setEventFileRegistry
(
    List<patchEventFile*>* eventFileRegistry,
    const word& dtFieldNameOverride
)
{
    eventFileRegistry_ = eventFileRegistry;
    dtFieldNameOverride_ = dtFieldNameOverride;
}

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
        //- Check if patchEventFile
        if (eventFileRegistry_)
        {
            eventFileRegistry_->append(&eventFile_);
        }
        else
        {
            FatalErrorIn("eventFlux.C")
                << "eventFlux BC is used with an incompatible solver"
                    << abort(FatalError);
        }

        //- reading patch event file, compute current value, store to old values
        eventFile_.read(eventFileName,true);
        eventFile_.updateIndex(this->db().time().startTime().value());
        eventFile_.storeOldValues();

        const word& dtFieldName = dtFieldNameOverride_.empty() ? iF.name() : dtFieldNameOverride_;
        eventFile_.setTimeScheme(dtFieldName, iF.mesh());

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
        valueEvent = eventFile_.dtValue(patchEventID_);
    }

    if ((mag(valueEvent + eventFluxValue_) > SMALL) && (mag(sum(phip_)) < VSMALL))
    {
        FatalErrorIn("eventFlux.C")
            << "non-zero fixed flux for C with zero flux field phi" << abort(FatalError);
    }

    //- Computing fixed value
    scalarField results(patch().patch().faceCentres().size());    
    results = (eventFluxValue_+valueEvent)/sum(phip_+VSMALL);
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
