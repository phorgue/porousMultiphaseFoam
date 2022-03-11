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

#include "XYfile.H"
#include "IFstream.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XYfile::XYfile
(
    const XYfile& XYtoCopy
)
:  
    name_(XYtoCopy.name()),
    mesh_(XYtoCopy.mesh()),
    x_(XYtoCopy.x()),
    y_(XYtoCopy.y()),
    values_(XYtoCopy.values()),
    mapping_(XYtoCopy.mapping())
{
}

Foam::XYfile::XYfile
(
    const word& fileName,
    const fvMesh& mesh,
    const label npoints
)
    :
    name_(fileName),
    mesh_(mesh),
    npoints_(npoints),
    mapping_(mesh.C().size())
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

    constructMapping();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XYfile::~XYfile()
{}

// * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

void Foam::XYfile::findClosestPoints(const point& location, labelList& id, scalarList& coeffs)
{
    id = -1;
    scalarList dist(npoints_, GREAT);
    forAll(x_, pointi)
    {
        scalar current_dist = Foam::sqrt(pow(x_[pointi]-location.x(),2)+pow(y_[pointi]-location.y(),2));

        //- finding point position in distance list
        label position = 0;
        forAll(dist, close_pointi)
        {
            if (current_dist > dist[close_pointi]) position++;
        }
        if (position < npoints_)
        {
            for(label iter=npoints_-1;iter>position;iter--)
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

    coeffs = 1/dist;
    scalar sumCoeffs = sum(coeffs);
    forAll(coeffs, pointi) coeffs[pointi] /= sumCoeffs;
}

void Foam::XYfile::constructMapping()
{
    word mappingFileName = name_ + ".map";
    OFstream mappingFile(mappingFileName);
    mappingFile << mesh_.C().size() << " " << npoints_ << endl;
    forAll(mesh_.C(), celli)
    {
        mapping_[celli].resize(npoints_);
        labelList id(npoints_);
        scalarList dist(npoints_);
        findClosestPoints(mesh_.C()[celli], id, dist);
        for(label i=0;i<npoints_;i++)
        {
            mapping_[celli][i].first() = id[i];
            mapping_[celli][i].second() = dist[i];
            mappingFile << mapping_[celli][i].first() << " " << mapping_[celli][i].second() << " ";
        }
        mappingFile << endl;
    };
}
void Foam::XYfile::readMapping()
{
}


// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //
Foam::scalar Foam::XYfile::interpolate(const point& location, const scalar& offset)
{
    labelList id(npoints_);
    scalarList coeffs(npoints_);
    findClosestPoints(location, id, coeffs);

    scalar interpolatedValue_ = offset;
    forAll(id,pointi) interpolatedValue_ += coeffs[pointi]*values_[id[pointi]];
    return interpolatedValue_;
}

void Foam::XYfile::mapField(volScalarField& field, const scalar& offset)
{
    forAll(field, celli)
    {
        field[celli] = offset;
        forAll(mapping_[celli], pointi)
        {
            field[celli] += values_[mapping_[celli][pointi].first()] * mapping_[celli][pointi].second();
        }
    }
}
