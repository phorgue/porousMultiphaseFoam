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

#include "patchEventFile.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::patchEventFile::patchEventFile()
    :
    eventFile(""),
    npatches_(0),
    patchNameList_()
{
}

Foam::patchEventFile::patchEventFile
(
    const patchEventFile& eventFileToCopy
)
    :
    eventFile(eventFileToCopy),
    npatches_(eventFileToCopy.npatches()),
    patchNameList_(eventFileToCopy.patchNameList())
{
}

Foam::patchEventFile::patchEventFile
(
    const word& fileName,
    bool display
)
    :
    eventFile(fileName)
{
    this->read(fileName,display);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchEventFile::~patchEventFile()
{}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //

void Foam::patchEventFile::read(const word& fileName, bool display)
{
    name_ = fileName;
    if (fileName.size() != 0)
    {
        //- properties of a MNT file
        string separator_ = " ";

        //- file name
        IFstream ifs(fileName);
        DynamicList<scalar> datesRead;
        DynamicList<word> nameRead;
        DynamicList<scalar> valueRead;
        if (display)
        {
            Info << nl << "Reading PatchEvent file '" << fileName << "' ...";
        }

        // read data
        while (ifs.good())
        {
            string line;
            ifs.getLine(line);

            if (line != "")
            {

                label n = 0;
                std::size_t pos = 0;
                DynamicList<string> split;

                while (pos != std::string::npos)
                {
                    std::size_t nPos = line.find(separator_, pos);
                    if (nPos == std::string::npos)
                    {
                        if (line.substr(pos).size() != 0)
                        {
                            split.append(line.substr(pos));
                            n++;
                        }
                        pos = nPos;
                    }
                    else
                    {
                        if (nPos - pos != 0)
                        {
                            split.append(line.substr(pos, nPos - pos));
                            n++;
                        }
                        pos = nPos + 1;
                    }
                }

                if (split.size() <= 1)
                {
                    break;
                }            
                else if (split[0] == "date")
                {
                    scalar newDate = readScalar(IStringStream(split[1])());
                    datesRead.append(newDate);
                }
                else
                {
                    if (n != 2)
                    {
                        FatalErrorIn("patchEventFile.C")
                            << "wrong number of elements in patchEvent file :" << fileName
                                << nl << " found " << split.size() << " elements instead of 2 "
                                << nl << "List of read elements : " << split
                                << abort(FatalError);
                    }

                    word patchName = IStringStream(split[0])();
                    nameRead.append(patchName);
                    scalar value = readScalar(IStringStream(split[1])());
                    valueRead.append(value);
                }
            }
        }

        ndates_ = datesRead.size();
        npatches_ = nameRead.size()/ndates_;

        if (display)
        {
            Info << "OK!"
                << nl << "{"
                << nl << "  number of dates   = " << ndates_
                << nl << "  number of patches = " << npatches_
                << nl << "  number of datas   = " << valueRead.size()
                << nl << "}" << endl;
        }

        //- Storing dates
        dates_.resize(ndates_);
        forAll(datesRead,datei)
        {
            dates_[datei] = datesRead[datei];
        }

        //- Storing patch name
        patchNameList_.resize(npatches_);
        forAll(patchNameList_,patchi)
        {
            patchNameList_[patchi] = nameRead[patchi];
        }
 
        //- Storing infiltration datas
        datas_.setSize(ndates_,npatches_);
        label iter = 0;
        forAll(datesRead,datei)
        {
            for(label patchi=0;patchi<npatches_;patchi++)
            {           
                datas_[datei][patchi] = valueRead[iter];
                iter++;
            }
        }

        currentValues_.setSize(npatches_);
        oldValues_.setSize(npatches_);
    }
}
