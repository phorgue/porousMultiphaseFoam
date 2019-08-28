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

#include "eventInfiltration.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eventInfiltration::
eventInfiltration
(
    const fvPatch& h,
    const DimensionedField<vector, volMesh>& iF
)
    :
    fixedValueFvPatchVectorField(h, iF),
    fixedInfiltrationValue_(0.),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{}


Foam::eventInfiltration::
eventInfiltration
(
    const fvPatch& h,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedValueFvPatchVectorField(h, iF, dict, false),
    fixedInfiltrationValue_(dict.lookupOrDefault<scalar>("constantValue",0.)),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{
    word eventFileName = dict.lookupOrDefault<word>("eventFile","");
    Info << nl << "eventFileName " << eventFileName << endl;
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
            FatalErrorIn("eventInfiltration.C") << " patch '" << patch().name() << "' not found in event file : " << eventFile_.name() << abort(FatalError);
        }
    }
    else
    {
        Info << "eventInfiltration boundary condition without event file" << endl;
    }

    //- Read if backward time scheme is used
    if (word(internalField().mesh().ddtScheme("source")) == "backward")
    {
        isBackwardScheme_ = true;
    }

}


Foam::eventInfiltration::
eventInfiltration
(
    const eventInfiltration& ptf,
    const fvPatch& h,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, h, iF, mapper),
    fixedInfiltrationValue_(ptf.fixedInfiltrationValue_),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{}


Foam::eventInfiltration::
eventInfiltration
(
    const eventInfiltration& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    fixedInfiltrationValue_(ptf.fixedInfiltrationValue_),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{}


Foam::eventInfiltration::
eventInfiltration
(
    const eventInfiltration& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    fixedInfiltrationValue_(ptf.fixedInfiltrationValue_),
    isBackwardScheme_(false),
    patchEventID_(-1),
    eventFile_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eventInfiltration::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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
    
    vectorField updatedFixedValue = (fixedInfiltrationValue_+valueEvent)*patch().nf();
    operator== (updatedFixedValue);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::eventInfiltration::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchVectorField,
    eventInfiltration
);
}

// ************************************************************************* //
