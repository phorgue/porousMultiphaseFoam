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

#include "outputEventFile.H"
#include "IFstream.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::outputEventFile> Foam::outputEventFile::New(Foam::Time& runTime)
{
    const bool isPresent = runTime.controlDict().found("outputEventFile");
    const bool CSVoutput = runTime.controlDict().found("CSVoutput");
    word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("outputEventFile","");
    return autoPtr<Foam::outputEventFile>(new outputEventFile(outputEventFileName, isPresent, CSVoutput, runTime));
}

Foam::outputEventFile::outputEventFile
(
    const word& fileName,
    const bool& isPresent,
    const bool& CSVoutput,
    Time& runTime
)
    :
    eventFile(fileName),
    isPresent_(isPresent),
    CSVoutput_(CSVoutput),
    runTime_(runTime),
    zscale_(1),
    fieldsToWrite_(0),
    coeff1Fields_(0),
    coeff2Fields_(0),
    phiFields_(0),
    massBalance_(0)
{
    if (isPresent_)
    {
        checkControlDict();
        if (fileName.size() != 0)
        {
            //- properties of a DEM file
            string separator_ = " ";

            //- file name
            IFstream ifs(fileName);
            DynamicList<scalar> datesRead;

            Info << nl << "Reading output event file '" << fileName << "' ...";
            // read data
            while (ifs.good()) {
                string line;
                ifs.getLine(line);

                if (line != "") {

                    label n = 0;
                    std::size_t pos = 0;
                    DynamicList<string> split;

                    while (pos != std::string::npos) {
                        std::size_t nPos = line.find(separator_, pos);
                        if (nPos == std::string::npos) {
                            if (line.substr(pos).size() != 0) {
                                split.append(line.substr(pos));
                                n++;
                            }
                            pos = nPos;
                        } else {
                            if (nPos - pos != 0) {
                                split.append(line.substr(pos, nPos - pos));
                                n++;
                            }
                            pos = nPos + 1;
                        }
                    }
                    if (n == 1) {
                        scalar newDate = readScalar(IStringStream(split[0])());
                        datesRead.append(newDate);
                    } else {
                        FatalErrorIn("outputEventFile.C")
                                << "wrong number of elements in event file :" << fileName
                                << nl << " found " << split.size() << " elements instead of 1 "
                                << nl << "List of read elements : " << split
                                << abort(FatalError);
                    }
                }
            }

            Info << "OK" << endl;
            ndates_ = datesRead.size();

            if (ndates_ == 0) {
                FatalErrorIn("outputEventFile.C")
                        << "Zero date in file : " << fileName
                        << nl << "File is empty or not present"
                        << abort(FatalError);
            }

            //- Storing dates
            dates_.resize(ndates_);
            forAll(datesRead, datei) {
                dates_[datei] = datesRead[datei];
            }

        }
        updateIndex(runTime_.startTime().value());
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::outputEventFile::~outputEventFile()
{}
// * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::outputEventFile::computeInterpolationFactor()
{
    scalar ifactor = (runTime_.timeOutputValue()-currentEventEndTime())/runTime_.deltaTValue();
    if (ifactor < 0 || ifactor > 1)
    {
        FatalErrorIn("outputEventFile.C")
            << " Unconsistent value for time interpolation = " << ifactor
            << " current time is " << runTime_.timeOutputValue()
            << " and event end time is " << currentEventEndTime() << abort(FatalError);
    }
    return ifactor;
}

// * * * * * * * * * * * * * * * * Members * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::outputEventFile::timeInterpolate
(
    const scalar& prev,
    const scalar& current
)
{
    //- compute interpolation factor
    scalar interpolateFactor = computeInterpolationFactor();
    return  interpolateFactor*current+(1.0-interpolateFactor)*prev;
}


template<class Type, template<class> class PatchField, class TypeMesh>
Foam::GeometricField<Type, PatchField, TypeMesh> Foam::outputEventFile::timeInterpolate
(
    const GeometricField<Type, PatchField, TypeMesh>& vfield,
    bool writeField
)
{
    //- compute interpolation factor
    scalar interpolateFactor = computeInterpolationFactor();

    //- update time
    scalar timeOutputBackup = runTime_.timeOutputValue();
    runTime_.setTime(currentEventEndTime(), runTime_.timeIndex());

    GeometricField<Type, PatchField, TypeMesh> ifield
        (
            IOobject
            (
                vfield.name(),
                runTime_.timeName(),
                vfield.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vfield
        );
    ifield = interpolateFactor*vfield+(1.0-interpolateFactor)*vfield.oldTime();
    if (writeField) ifield.write();

    runTime_.setTime(timeOutputBackup,runTime_.timeIndex());
    return ifield;
}

void Foam::outputEventFile::checkControlDict() const
{
    const dictionary& cDict = runTime_.controlDict();
    scalar endTimeValue(cDict.get<scalar>("endTime"));
    scalar writeIntervalValue(cDict.getOrDefault<scalar>("writeInterval",GREAT));
    if (cDict.found("writeFrequency"))
    {
        FatalErrorIn("outputEventFile.C") << "Check controlDict: writeFrequency cannot be used when outputEventFile is specified"<< exit(FatalError);
    }
    if (endTimeValue > writeIntervalValue)
    {
        FatalErrorIn("outputEventFile.C") << "Check controlDict: writeInterval should be larger or equal to endTime when outputEventFile is specified"<< exit(FatalError);
    }

}
void Foam::outputEventFile::addField(
    const volScalarField& field,
    const surfaceScalarField& phi,
    const word& type,
    bool massBalance,
    const label nCoef
) {
    if (nCoef<2) coeff2Fields_.append(nullptr);
    if (nCoef<1) coeff1Fields_.append(nullptr);
    fieldsToWrite_.append(&field);
    phiFields_.append(&phi);
    massBalance_.append(massBalance);
    if (CSVoutput_ && massBalance)
    {
        CSVoutputFiles_.append(new OFstream(field.name() + "massBalance.csv"));
        OFstream& massBalanceCSV = CSVoutputFiles_[CSVoutputFiles_.size()-1];
        massBalanceCSV << "#Time Total(" << type << ")";
        const fvMesh& mesh = field.mesh();
        forAll(mesh.boundaryMesh(),patchi)
        {
            if (mesh.boundaryMesh()[patchi].type() == "patch")
            {
                massBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
            }
        }
        massBalanceCSV << endl;
    }
    else
    {
        CSVoutputFiles_.append(nullptr);
    }
}

void Foam::outputEventFile::addField(
    const volScalarField& field,
    const surfaceScalarField& phi,
    const volScalarField& coef1,
    const word& type,
    bool massBalance
) {
    coeff1Fields_.append(&coef1);
    addField(field, phi, type, massBalance, 1);
}

void Foam::outputEventFile::addField(
    const volScalarField& field,
    const surfaceScalarField& phi,
    const volScalarField& coef1,
    const volScalarField& coef2,
    const word& type,
    bool massBalance
) {
    coeff1Fields_.append(&coef1);
    coeff2Fields_.append(&coef2);
    addField(field, phi, type, massBalance, 2);
}

void Foam::outputEventFile::addField(
    const volScalarField& field,
    const surfaceScalarField& phi,
    const volScalarField& coef1,
    const volScalarField& coef2,
    const volScalarField& coef3,
    const word& type,
    bool massBalance
) {
    coeff1Fields_.append(&coef1);
    coeff2Fields_.append(&coef2);
    coeff3Fields_.append(&coef3);
    addField(field, phi, type, massBalance, 3);
}

void Foam::outputEventFile::write() {
    if (isPresent_) {
        if (currentEventEndTime() <= runTime_.timeOutputValue()) {
            forAll(fieldsToWrite_, fieldi) {
                const fvMesh& mesh = fieldsToWrite_[fieldi].mesh();
                volScalarField fInter = timeInterpolate(fieldsToWrite_[fieldi], true);
                surfaceScalarField phiInter =  timeInterpolate(phiFields_[fieldi], true);

                if (CSVoutput_ && massBalance_[fieldi]) {
                    OFstream& massBalanceCSV = CSVoutputFiles_[fieldi];
                    if (coeff1Fields_.test(fieldi)) {
                        volScalarField cInter1 = timeInterpolate(coeff1Fields_[fieldi], false);
                        if (coeff2Fields_.test(fieldi)) {
                            volScalarField cInter2 = timeInterpolate(coeff1Fields_[fieldi], false);
                            if (coeff3Fields_.test(fieldi)) {
                                volScalarField cInter3 = timeInterpolate(coeff1Fields_[fieldi], false);
                                massBalanceCSV << currentEventEndTime() << " " << fvc::domainIntegrate(fInter * cInter1 *  cInter2 * cInter3).value()/zscale_;
                            }
                            else {
                                massBalanceCSV << currentEventEndTime() << " " << fvc::domainIntegrate(fInter * cInter1 *  cInter2).value()/zscale_;
                            }
                        }
                        else {
                            massBalanceCSV << currentEventEndTime() << " " << fvc::domainIntegrate(fInter * cInter1).value()/zscale_;
                        }
                    }
                    else {
                        massBalanceCSV << currentEventEndTime() << " " << fvc::domainIntegrate(fInter).value()/zscale_;
                    }
                    forAll(mesh.boundaryMesh(), patchi) {
                        if (mesh.boundaryMesh()[patchi].type() == "patch") {
                            massBalanceCSV << " "
                                           << gSum(phiInter.boundaryField()[patchi] * fInter.boundaryField()[patchi]);
                        }
                    }
                    massBalanceCSV << endl;
                }
            }
            updateIndex(runTime_.timeOutputValue());

        }
    }
    else {
        runTime_.write();
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
