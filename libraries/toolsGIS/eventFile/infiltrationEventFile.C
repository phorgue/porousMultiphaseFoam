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

#include "infiltrationEventFile.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::infiltrationEventFile::infiltrationEventFile
(
    const infiltrationEventFile& eventFileToCopy
)
    :
    eventFile(eventFileToCopy),
    uniform_(eventFileToCopy.uniform_)
{
}

Foam::infiltrationEventFile::infiltrationEventFile
(
    const word& fileName
)
    :
    eventFile(fileName),
    uniform_(true)
{
    if (fileName.size() != 0)
    {
        //- Infiltration field size (1 = uniform infiltration)
        label fieldSize = 1;

        //- properties of a DEM file
        string separator_ = " ";

        //- file name
        IFstream ifs(fileName);
        DynamicList<scalar> datesRead;
        DynamicList<DynamicList<scalar> > valueRead;

       // read data
        Info << nl << "Reading Event file '" << fileName << "' ...";
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

                if (split.size() < 1)
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
                    DynamicList<scalar> tmpValueRead;
                    for(label iter=0;iter<split.size();iter++)
                    {
                        scalar newValue =  readScalar(IStringStream(split[iter])());
                        tmpValueRead.append(newValue);
                    }
                    if (tmpValueRead.size() > 1)
                    {
                        if (fieldSize == 1) fieldSize = tmpValueRead.size();
                        if (fieldSize != tmpValueRead.size())
                        {
                            FatalErrorIn("infiltrationEventFile.C")
                                << "wrong number of elements in event file :" << fileName
                                    << nl << " found " << split.size() << " elements instead of 1 or "
                                    << fieldSize << " (size of the first non-uniform infiltration data)"
                                    << abort(FatalError);
                        }
                    }
                    valueRead.append(tmpValueRead);
                }

            }
        }
        ndates_ = datesRead.size();
    
        Info << "OK!"
            << nl << "{"
            << nl << "  number of dates      = " << ndates_;
        if (fieldSize == 1)
        {
            Info << nl << "  type of infiltration = uniform";
        }
        else
        {
            Info << nl << "  type of infiltration = nonuniform";
        }
        Info << nl << "}" << endl;

        //- Storing dates
        dates_.resize(ndates_);
        forAll(datesRead,datei)
        {
            dates_[datei] = datesRead[datei];
        }

        //- Storing infiltration datas
        datas_.setSize(ndates_,fieldSize);
        forAll(datesRead,datei)
        {
            scalarList currentData = valueRead[datei];
            if (currentData.size() == 1)
            {
                for(label celli=0;celli<fieldSize;celli++)
                {
                    datas_[datei][celli] = currentData[0];
                }
            }
            else
            {
                for(label celli=0;celli<fieldSize;celli++)
                {
                    datas_[datei][celli] = currentData[celli];
                }
            }
        }

        if (fieldSize > 1) uniform_ = false;
        currentValues_.setSize(fieldSize,0);
        oldValues_.setSize(fieldSize,0);
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::infiltrationEventFile::~infiltrationEventFile()
{}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //
