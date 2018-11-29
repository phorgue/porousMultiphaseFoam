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

#include "MNTfile.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MNTfile::MNTfile
(
    const MNTfile& MNTtoCopy
)
:  
    name_(MNTtoCopy.name()),
    x0_(MNTtoCopy.x0()),
    y0_(MNTtoCopy.y0()),
    dx_(MNTtoCopy.dx()),
    dy_(MNTtoCopy.dy()),
    zvalues_(MNTtoCopy.z())
{
}

Foam::MNTfile::MNTfile
(
    const word& fileName
)
    :
    name_(fileName)
{
    //- properties of a MNT file
    string separator_ = " ";
    label nEntries = 3;

    //- file name
    IFstream ifs(fileName);
    DynamicList<point> mnt_data; 

    // read data
    while (ifs.good())
    {
        string line;
        ifs.getLine(line);

        label n = 0;
        std::size_t pos = 1;
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

        if (n != 3)
        {
            FatalErrorIn("MNTfile.C")
                << "wrong number of elements in MNT file :" << fileName
                    << nl << " found " << split.size() << " elements instead of 3 atline " << mnt_data.size()+1
                    << nl << "List of read elements : " << split
                    << abort(FatalError);
        }

        scalar x = readScalar(IStringStream(split[0])());
        scalar y = readScalar(IStringStream(split[1])());
        scalar z = readScalar(IStringStream(split[2])());
        mnt_data.append(point(x,y,z));

    }

    Info << nl << "MNT file : " << fileName << " read and found " << mnt_data.size() << " points" << endl;
    
    x0_ = mnt_data[0][0];
    y0_ = mnt_data[0][1];
    dx_ = mnt_data[1][0] - mnt_data[0][0];
    dy_ = 0;
    if (dx_ <= 0) FatalErrorIn("MNTfile.C") << "MNT file incorrectly sorted" << abort(FatalError);
    zvalues_ = scalarList(mnt_data.size());
    zvalues_[0] = mnt_data[0][2];

    for (label i=1;i<mnt_data.size();i++)
    {
        zvalues_[i] = mnt_data[i][2];
        if (dy_ == 0)
        {
            if (mnt_data[i][1] > y0_)
            {
                dy_ = mnt_data[i][1] - y0_;
                nx_ = i+1;
            }
        }
    }
    ny_ = mnt_data.size()/nx_;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MNTfile::~MNTfile()
{}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //
Foam::scalar Foam::MNTfile::interpolate(const point& location)
{
    label idx_ = floor((location.x()-x0_)/dx_);
    label idy_ = floor((location.y()-y0_)/dy_);
    scalar fracx_ = (location.x() - idx_*dx_ - x0_)/dx_;
    scalar fracy_ = (location.y() - idy_*dy_ - y0_)/dy_;
    scalar tmp1_ = (1.0-fracx_) * zvalues_[idx_+nx_*idy_]+ fracx_ * zvalues_[idx_+1+nx_*idy_];
    scalar tmp2_ = (1.0-fracx_) * zvalues_[idx_+nx_*(idy_+1)]+ fracx_ * zvalues_[idx_+1+nx_*(idy_+1)];
    scalar interpolatedValue_ = (1.0-fracy_) * tmp1_ + fracy_ * tmp2_;
    return interpolatedValue_;
}
