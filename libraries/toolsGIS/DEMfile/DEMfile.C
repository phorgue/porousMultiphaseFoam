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

#include "DEMfile.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DEMfile::DEMfile
(
    const DEMfile& DEMtoCopy
)
:  
    name_(DEMtoCopy.name()),
    x0_(DEMtoCopy.x0()),
    y0_(DEMtoCopy.y0()),
    dx_(DEMtoCopy.dx()),
    dy_(DEMtoCopy.dy()),
    zvalues_(DEMtoCopy.z())
{
}

Foam::DEMfile::DEMfile
(
    const word& fileName
)
    :
    name_(fileName)
{
    //- properties of a DEM file
    string separator_ = " ";
    label nEntries = 3;

    //- file name
    IFstream ifs(fileName);
    DynamicList<point> mnt_data; 

    Info << nl << "Reading DEM file '" << fileName << "' ...";
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

            if (n != 3)
            {
                FatalErrorIn("DEMfile.C")
                    << "wrong number of elements in DEM file :" << fileName
                        << nl << " found " << split.size() << " elements instead of 3 at line " << mnt_data.size()+1
                        << nl << "List of read elements : " << split
                        << abort(FatalError);
            }

            scalar x = readScalar(IStringStream(split[0])());
            scalar y = readScalar(IStringStream(split[1])());
            scalar z = readScalar(IStringStream(split[2])());
            mnt_data.append(point(x,y,z));
        }

    }

    Info << "OK!"
        << nl << "{"
        << nl << "  number of points = " << mnt_data.size() << " ";

    x0_ = mnt_data[0][0];
    y0_ = mnt_data[0][1];
    dx_ = mnt_data[1][0] - mnt_data[0][0];
    dy_ = 0;
    if (dx_ <= 0) FatalErrorIn("DEMfile.C") << "DEM file incorrectly sorted (x1 = x2)" << abort(FatalError);
    zvalues_ = scalarList(mnt_data.size());
    zvalues_[0] = mnt_data[0][2];

    for (label i=1;i<mnt_data.size();i++)
    {
        zvalues_[i] = mnt_data[i][2];
        if (dy_ == 0)
        {
            if (mnt_data[i][1] != y0_)
            {
                dy_ = mnt_data[i][1] - y0_;
                nx_ = i;
            }
        }
    }
    ny_ = mag(mnt_data.size()/nx_);

    Info << nl << "  grid = " << nx_ << " x " << ny_
        << nl << "  dx = " << dx_ << " and dy = " << dy_
        << nl << "  startPoint (" << x0_ << "," << y0_ << ")"
        << nl << "  endPoint   (" << x0_+nx_*dx_ << "," << y0_+ny_*dy_ << ")"
        << nl << "}" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DEMfile::~DEMfile()
{}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //
Foam::scalar Foam::DEMfile::interpolate(const point& location)
{
    label idx_ = floor((location.x()-x0_)/dx_);
    label idy_ = floor((location.y()-y0_)/dy_);
    if ((idx_ < 0) || (idx_ >= nx_) || (idy_ < 0) || (idy_ >= ny_))
    {
        FatalErrorIn("Foam::scalar Foam::DEMfile::interpolate(const point& location)") << "location "
            << location << " out of DEM bounds" << abort(FatalError);
    }
    scalar fracx_ = (location.x() - idx_*dx_ - x0_)/dx_;
    scalar fracy_ = (location.y() - idy_*dy_ - y0_)/dy_;
    scalar tmp1_ = (1.0-fracx_) * zvalues_[idx_+nx_*idy_]+ fracx_ * zvalues_[idx_+1+nx_*idy_];
    scalar tmp2_ = (1.0-fracx_) * zvalues_[idx_+nx_*(idy_+1)]+ fracx_ * zvalues_[idx_+1+nx_*(idy_+1)];
    scalar interpolatedValue_ = (1.0-fracy_) * tmp1_ + fracy_ * tmp2_;
    return interpolatedValue_;
}
