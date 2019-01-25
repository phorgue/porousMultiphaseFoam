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

Class
    Foam::MNTfile

Description
    MNT file handler/reader with bilinear interpolation for the
    porousMultiphaseFoam toolbox

SourceFiles
    MNTfile.C

\*---------------------------------------------------------------------------*/

#ifndef MNTfile_H
#define MNTfile_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class fluidMNTfile Declaration
    \*---------------------------------------------------------------------------*/

class MNTfile
{

private:

    word name_;

    scalar x0_,y0_;
    scalar dx_,dy_;
    label nx_,ny_;
    scalarList zvalues_;

public:

    // Constructors

    //- Copy constructor
    MNTfile
    (
        const MNTfile& fileToCopy
    );

    //- Construct from Istream
    MNTfile
    (
        const word& fileName
    );
    
    // Selectors

    //- Destructor
    virtual ~MNTfile();

    // Member Functions
    
    //- Access function
    const word& name() const {return name_;}
    const scalar& x0() const {return x0_;}
    const scalar& y0() const {return y0_;}
    const scalar& dx() const {return dx_;}
    const scalar& dy() const {return dy_;}
    const scalarList& z() const {return zvalues_;}

    //- Bilinear interpolation using given point
    scalar interpolate(const point& position);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
