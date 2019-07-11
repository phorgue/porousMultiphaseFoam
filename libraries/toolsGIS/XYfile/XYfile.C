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

#include "XYfile.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XYfile::XYfile
(
    const XYfile& XYtoCopy
)
:  
    name_(XYtoCopy.name()),
    x_(XYtoCopy.x()),
    y_(XYtoCopy.y()),
    values_(XYtoCopy.values())
{
}

Foam::XYfile::XYfile
(
    const word& fileName
)
    :
    name_(fileName)
{
    //- properties of a XY file
    string separator_ = " ";
    label nEntries = 3;

    //- file name
    IFstream ifs(fileName);
    DynamicList<scalar> xread; 
    DynamicList<scalar> yread;
    DynamicList<scalar> valuesread;

    Info << nl << "Reading XY file '" << fileName << "' ...";
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
        
            while ((pos != std::string::npos) && (n <= nEntries))
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

            if (n != nEntries)
            {
                FatalErrorIn("XYfile.C")
                    << "wrong number of elements in XY file :" << fileName
                    //        << nl << " found " << split.size() << " elements instead of 3 at line " << mnt_data.size()+1
                        << nl << "List of read elements : " << split
                        << abort(FatalError);
            }

            xread.append(readScalar(IStringStream(split[0])()));
            yread.append(readScalar(IStringStream(split[1])()));
            valuesread.append(readScalar(IStringStream(split[2])()));
            
        }
    }

    //- storing data read
    x_.resize(xread.size());
    x_=xread;
    y_.resize(yread.size());
    y_=yread;
    values_.resize(valuesread.size());
    values_=valuesread;

    //- display general informations
    Info << "OK!"
        << nl << "{"
        << nl << "  number of points = " << x_.size() << " "
        << nl << "  startPoint (" << min(x_) << "," << min(y_) << ")"
        << nl << "  endPoint   (" << max(x_) << "," << max(y_) << ")"
        << nl << "}" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XYfile::~XYfile()
{}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //
Foam::scalar Foam::XYfile::interpolate(const point& location)
{
    label id1=-1;
    label id2=-1;
    label id3=-1;     
    scalar dist1 = GREAT;
    scalar dist2 = GREAT;
    scalar dist3 = GREAT;

    for(label i=0;i<x_.size();i++)
    {
        scalar dist = Foam::sqrt(pow(x_[i]-location.x(),2)+pow(y_[i]-location.y(),2));
        if ( dist < dist1)
        {
             id3 = id2;
            dist3 = dist2;
            id2 = id1;
            dist2 = dist1;
            id1 = i;
            dist1 = dist;
        }
        else if ( dist < dist2)
        {   
            id3 = id2;
            dist3 = dist2;
            id2 = i;
            dist2 = dist;                
        }
        else if ( dist < dist3)
        {
                 id3 = i;
            dist3 = dist;                
        }
    }

        
    if ( (id1 == -1) || (id2 == -1) || (id3 == -1))
    {
        Info << nl << "Error : three point are not found for interpolation"
            << nl << id1 << " / " << id2 << " / " << id3 << endl;
    }

    scalar interpolatedValue_ =  (dist1*values_[id1] + dist2*values_[id2] + dist3*values_[id3] ) / (dist1+dist2+dist3);
    return interpolatedValue_;
}
