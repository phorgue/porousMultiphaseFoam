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

#include "eventInfiltration.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"


Foam::List<Foam::patchEventFile*>* Foam::eventInfiltration::eventFileRegistry_ = nullptr;
Foam::word Foam::eventInfiltration::dtFieldNameOverride_ = "";

void Foam::eventInfiltration::setEventFileRegistry
(
    List<patchEventFile*>* eventFileRegistry,
    const word& dtFieldNameOverride
)
{
    eventFileRegistry_ = eventFileRegistry;
    dtFieldNameOverride_ = dtFieldNameOverride;
}

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
    fixedInfiltrationValue_(dict.getOrDefault<scalar>("constantValue",0.)),
    patchEventID_(-1),
    eventFile_()
{
    word eventFileName = dict.getOrDefault<word>("eventFile","");
    Info << nl << "eventFileName " << eventFileName << endl;

    if (eventFileName != "")
    {
        //- Check if patchEventFile
        if (eventFileRegistry_)
        {
            eventFileRegistry_->append(&eventFile_);
        }
        else
        {
            FatalErrorIn("eventInfiltration.C")
                << "eventInfiltration BC is used with an incompatible solver"
                    << abort(FatalError);
        }

        //- reading patch event file, compute current value, store to old values
        eventFile_.read(eventFileName,true);
        eventFile_.updateIndex(this->db().time().startTime().value());
        eventFile_.storeOldValues();

        const word& dtFieldName = dtFieldNameOverride_.empty() ? iF.name() : dtFieldNameOverride_;
        eventFile_.setTimeScheme(dtFieldName, iF.mesh());

        //- Reading patch event file and adding intermediate time step
        scalar eventTimeStep = this->db().time().controlDict().getOrDefault<scalar>("eventTimeStep",0);
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
        valueEvent = eventFile_.dtValue(patchEventID_);
    }

    //- Computing fixed value
    
    vectorField updatedFixedValue((fixedInfiltrationValue_+valueEvent)*patch().nf());
    operator== (updatedFixedValue);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::eventInfiltration::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    this->writeEntry("value", os);
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
