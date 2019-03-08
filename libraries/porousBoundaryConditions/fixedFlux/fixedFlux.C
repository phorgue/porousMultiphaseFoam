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
    word eventFileName = db().lookupObject<dictionary>("transportProperties").lookupOrDefault<word>("eventFilePatchMassFlowRate","");
    
    if (eventFileName != "")
    {
        //- Reading patch event file and adding intermediate time step
        patchEventFile eventFlux(eventFileName);
        scalar eventTimeStep = this->db().time().controlDict().lookupOrDefault<scalar>("eventTimeStep",0);
        if (eventTimeStep > 0)
        {
            eventFlux.addIntermediateTimeSteps(eventTimeStep);
        }

        //- finding patch name
        label patchID = -1;
        forAll(eventFlux.patchNameList(),patchi)
        {
            if (patch().name() == eventFlux.patchNameList()[patchi]) patchID = patchi;
        }
        if (patchID == -1)
        {
            FatalErrorIn("fixedFlux.C") << " patch '" << patch().name() << "' not found in event file : " << eventFlux.name() << abort(FatalError);
        }

        //- storing event data
        eventData_.setSize(eventFlux.ndates());
        forAll(eventData_,eventi)
        {
            eventData_[eventi] = Tuple2<scalar, scalar>(eventFlux.dates()[eventi],eventFlux.datas()[eventi][patchID]);
        }

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
    if (iterEvent_ > -1 )
    {
        scalar currentTime = this->db().time().value();
        if (currentTime > eventData_[iterEvent_+1].first())
        {
            iterEvent_++;
            scalar deltaT =this->db().time().deltaT().value();
            scalar deltaT0 =this->db().time().deltaT0().value();
            // Coefficient for t-3/2 (between times 0 and 00)
            scalar coefft0_00 = deltaT/(deltaT + deltaT0);

            // Coefficient for t-1/2 (between times n and 0)
            scalar coefftn_0 = 1 + coefft0_00;
            valueEvent_ = coefftn_0*eventData_[iterEvent_].second() - coefft0_00*eventData_[iterEvent_-1].second();
        }
        else
        {
        valueEvent_ = eventData_[iterEvent_].second();
        }
    }

    //- Computing fixed value
    scalarField results(patch().patch().faceCentres().size());    
    results = (fixedFluxValue_+valueEvent_)/sum(phip_);
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
