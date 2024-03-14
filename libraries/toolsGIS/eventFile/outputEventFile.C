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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::outputEventFile> Foam::outputEventFile::New(Foam::Time& runTime, const fvMesh& mesh, const scalar zscale)
{
    const bool isPresent = runTime.controlDict().found("outputEventFile");
    const bool CSVoutput = runTime.controlDict().lookupOrDefault("CSVoutput", false);
    word outputEventFileName = runTime.controlDict().lookupOrDefault<word>("outputEventFile","");
    return autoPtr<Foam::outputEventFile>(new outputEventFile(outputEventFileName, isPresent, CSVoutput, runTime, mesh, zscale));
}

Foam::outputEventFile::outputEventFile
(
    const word& fileName,
    const bool& isPresent,
    const bool& CSVoutput,
    Time& runTime,
    const fvMesh& mesh,
    const scalar zscale
)
    :
    eventFile(fileName),
    isPresent_(isPresent),
    CSVoutput_(CSVoutput),
    runTime_(runTime),
    zscale_(zscale),
    one_(
        IOobject
        (
            "one",
            runTime_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    ),
    outputFields_()
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

Foam::outputEventFile::~outputEventFile() = default;

// * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * * * * //

void Foam::outputEventFile::updateInterpolationFactor()
{
    ifactor_ = (runTime_.timeOutputValue()-currentEventEndTime())/runTime_.deltaTValue();
    if (ifactor_ < 0 || ifactor_ > 1)
    {
        FatalErrorIn("outputEventFile.C")
            << " Unconsistent value for time interpolation = " << ifactor_
            << " current time is " << runTime_.timeOutputValue()
            << " and event end time is " << currentEventEndTime() << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Members * * * * * * * * * * * * * * * * * //

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
    const volScalarField& coef1,
    const volScalarField& coef2,
    const volScalarField& coef3,
    word fileName,
    bool saturation
) {

    outputFields_.append(new outputField(
        field,
        phi,
        coef1,
        coef2,
        coef3,
        saturation,
        CSVoutput_,
        zscale_,
        fileName
        )
    );
}

void Foam::outputEventFile::addField(
    const volScalarField& field,
    const surfaceScalarField& phi,
    const volScalarField& coef1,
    const volScalarField& coef2,
    word fileName,
    bool saturation
) {
    addField(field, phi, coef1, coef2, one_, fileName, saturation);
}

void Foam::outputEventFile::addField(
        const volScalarField& field,
        const surfaceScalarField& phi,
        const volScalarField& coef1,
        word fileName,
        bool saturation
) {
    addField(field, phi, coef1, one_, one_, fileName, saturation);
}

void Foam::outputEventFile::addField(
        const volScalarField& field,
        const surfaceScalarField& phi,
        word fileName,
        bool saturation
) {
    addField(field, phi, one_, one_, one_, fileName, saturation);
}

void Foam::outputEventFile::write() {
    ifactor_ = 1;
    if (isPresent_) {
        if (currentEventEndTime() <= runTime_.timeOutputValue()) {
            updateInterpolationFactor();
            word timeEvent = Foam::name(currentEventEndTime());
            scalar timeBackup = runTime_.timeOutputValue();
            runTime_.setTime(currentEventEndTime(), runTime_.timeIndex());
            forAll(outputFields_, fieldi) {
                outputFields_[fieldi].write(timeEvent, ifactor_);
            }
            runTime_.setTime(timeBackup, runTime_.timeIndex());
            updateIndex(runTime_.timeOutputValue());
        }
    }
    else {
        runTime_.write();
        forAll(outputFields_, fieldi) {
            outputFields_[fieldi].write(runTime_.timeName());
        }
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
