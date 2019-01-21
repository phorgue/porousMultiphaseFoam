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

#include "eventFile.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eventFile::eventFile
(
    const eventFile& MNTtoCopy
)
    :  
    name_(MNTtoCopy.name()),
    ndates_(MNTtoCopy.ndates()),
    ncoordinates_(MNTtoCopy.ncoordinates()),
    dates_(MNTtoCopy.dates()),
    coordinates_(MNTtoCopy.coordinates()),
    datas_(MNTtoCopy.datas())
{
}

Foam::eventFile::eventFile
(
    const word& fileName
)
    :
    name_(fileName)
{
    if (fileName.size() != 0)
    {
        //- properties of a MNT file
        string separator_ = " ";
        label nEntriesMax = 4;

        //- file name
        IFstream ifs(fileName);
        DynamicList<scalar> datesRead;
        DynamicList<point> coordinatesRead;
        DynamicList<scalar> valueRead;

        Info << nl << "Reading Event file '" << fileName << "' ...";
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
        
                while ((pos != std::string::npos) && (n <= nEntriesMax))
                {
                    std::size_t nPos = line.find(separator_, pos);
                    if (nPos == std::string::npos)
                    {
                        split.append(line.substr(pos));
                        pos = nPos;
                        n++;
                    }
                    else
                    {
                        split.append(line.substr(pos, nPos - pos));
                        pos = nPos + 1;
                        n++;
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
                    if (n != 4)
                    {
                        FatalErrorIn("eventFile.C")
                            << "wrong number of elements in event file :" << fileName
                                << nl << " found " << split.size() << " elements instead of 4 "
                                << nl << "List of read elements : " << split
                                << abort(FatalError);
                    }
       
                    scalar x = readScalar(IStringStream(split[0])());
                    scalar y = readScalar(IStringStream(split[1])());
                    scalar z = readScalar(IStringStream(split[2])());
                    coordinatesRead.append(point(x,y,z));
                    scalar value = readScalar(IStringStream(split[3])());
                    valueRead.append(value);
                }
            }
        }

        ndates_ = datesRead.size();
        ncoordinates_ = coordinatesRead.size()/ndates_;
    
        Info << "OK!"
            << nl << "{"
            << nl << "  number of dates       = " << ndates_
            << nl << "  number of coordinates = " << ncoordinates_
            << nl << "  number datas          = " << coordinatesRead.size()
            << nl << "}" << endl;

        //- Storing dates
        dates_ = datesRead;
 
        //- Storing coordinates
        coordinates_.resize(ncoordinates_);
        forAll(coordinates_,coordinatei)
        {
            coordinates_[coordinatei] = coordinatesRead[1+coordinatei];
        }
 
        //- Storing infiltration datas
        datas_.setSize(ndates_,ncoordinates_);
        label iter = 0;
        forAll(datesRead,datei)
        {
            for(label coordinatei=0;coordinatei<ncoordinates_;coordinatei++)
            {           
                datas_[datei][coordinatei] = valueRead[iter];
                iter++;
            }
        }
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eventFile::~eventFile()
{}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //