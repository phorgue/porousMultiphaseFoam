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
    Foam::eventFile

Description
    Event file handler/reader for the porousMultiphaseFoam toolbox which
    contain 'm' dates with 'n' infiltration points where infiltration is
    expressed in m3/s. The class reads the file :

    date time1
    x1 y1 z1 infiltration11
    x2 y2 z2 infiltration12
    ...
    xn yn zn infiltration1n
    date time2
    x1 y1 z1 infiltration21
    x2 y2 z2 infiltration22
    ...
    date timem
    ...
    xn yn zn infiltrationmn

SourceFiles
    eventFile.C

Author
    Pierre Horgue

\*---------------------------------------------------------------------------*/

#ifndef eventFile_H
#define eventFile_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class fluideventFile Declaration
    \*---------------------------------------------------------------------------*/

class eventFile
{

private:

    word name_;
    label ndates_;
    label ncoordinates_;
    scalarList dates_;
    List<point> coordinates_;
    RectangularMatrix<scalar> datas_;

public:

    // Constructors    
    //- Copy constructor
    eventFile
    (
        const eventFile& fileToCopy
    );

    //- Construct from Istream
    eventFile
    (
        const word& fileName
    );
    
    // Selectors

    //- Destructor
    virtual ~eventFile();

    // Member Functions
    
    //- Access function
    const word& name() const {return name_;}
    const label& ndates() const {return ndates_;}
    const label& ncoordinates() const {return ncoordinates_;}
    const scalarList& dates() const {return dates_;}
    const List<point>& coordinates() const {return coordinates_;}
    const RectangularMatrix<scalar>& datas() const {return datas_;}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
