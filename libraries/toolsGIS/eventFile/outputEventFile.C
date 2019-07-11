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

#include "outputEventFile.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outputEventFile::outputEventFile
(
    const outputEventFile& eventFileToCopy
)
    :
    eventFile(eventFileToCopy)
{
}

Foam::outputEventFile::outputEventFile
(
    const word& fileName
)
    :
    eventFile(fileName)
{
    if (fileName.size() != 0)
    {
        //- properties of a MNT file
        string separator_ = " ";

        //- file name
        IFstream ifs(fileName);
        DynamicList<scalar> datesRead;

        Info << nl << "Reading output event file '" << fileName << "' ...";
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
                if (n == 1)
                {
                    scalar newDate = readScalar(IStringStream(split[0])());
                    datesRead.append(newDate);
                }
                else
                {
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
    
        //- Storing dates
        dates_.resize(ndates_);
        forAll(datesRead,datei)
        {
            dates_[datei] = datesRead[datei];
        }
    
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::outputEventFile::~outputEventFile()
{}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //
