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

Class
    Foam::XYfile

Description
    XY file handler/reader with linear interpolation for the
    porousMultiphaseFoam toolbox
    The file should be in the form

    x1 y1 value1
    x2 y2 value2
    ...
    xN yN valueN

SourceFiles
    XYfile.C

\*---------------------------------------------------------------------------*/

#ifndef XYfile_H
#define XYfile_H

#include "fvMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class XYfile Declaration
    \*---------------------------------------------------------------------------*/

class XYfile
{

private:

    word name_;
    const fvMesh& mesh_;
    label npoints_;
    scalarList x_;
    scalarList y_;
    scalarList values_;
    List<List<Tuple2<label, scalar>>> mapping_ ;

    DynamicList<string> splitLine(const string&, const string&, const label&, const word&);
    void findClosestPoints(const point&, labelList&, scalarList&);
    void constructMapping();
    bool readMapping(IFstream&);

public:

    // Constructors

    //- Copy constructor
    XYfile
    (
        const XYfile& fileToCopy
    );

    //- Construct from Istream
    XYfile
    (
        const word& fileName,
        const fvMesh& mesh,
        const label npoints=3
    );
    
    // Selectors

    //- Destructor
    virtual ~XYfile();

    // Member Functions
    
    //- Access function
    const word& name() const {return name_;}
    const fvMesh& mesh() const {return mesh_;}
    const scalarList& x() const {return x_;}
    const scalarList& y() const {return y_;}
    const scalarList& values() const {return values_;}
    const  List<List<Tuple2<label, scalar>>>& mapping() const {return mapping_;}

    //- Linear interpolation using given point (based on the 3 closest points)
    scalar interpolate(const point& position, const scalar& offset=0);
    void mapField(volScalarField& field, const scalar& offset=0);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
