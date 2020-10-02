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
Foam::scalar Foam::XYfile::interpolate(const point& location, label npoints=3)
{
    labelList id(npoints);
    id = -1;
    scalarList dist(npoints);
    dist = GREAT;

    for(label pointi=0;pointi<x_.size();pointi++)
    {
        scalar current_dist = Foam::sqrt(pow(x_[pointi]-location.x(),2)+pow(y_[pointi]-location.y(),2));

        //- finding point position in distance list
        label position = 0;
        forAll(dist, close_pointi)
        {
            if (current_dist > dist[close_pointi]) position++;
        }
        if (position < npoints)
        {
            for(label iter=npoints-1;iter>position;iter--)
            {
                id[iter] = id[iter-1];
                dist[iter] = dist[iter-1];
            }
            id[position] = pointi;
            dist[position] = current_dist+SMALL;
        }
    }
        
    if ( min(id) < 0 )
    {
        FatalErrorIn("XYfile.C") << nl << "Error : somepoint are not found for interpolation"
            << nl << id << abort(FatalError);
    }

    scalar interpolatedValue_ = 0;
    scalarList coeffs_(npoints);
    coeffs_ = 1/dist;
    scalar total_coeffs_ = sum(coeffs_);

    forAll(id,pointi) interpolatedValue_ += coeffs_[pointi]*values_[id[pointi]] / total_coeffs_;
    return interpolatedValue_;
}
