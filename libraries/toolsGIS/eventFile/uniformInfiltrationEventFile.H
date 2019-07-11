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
    Foam::uniformInfiltrationEventFile

Description
    Event file handler/reader for the porousMultiphaseFoam toolbox which
    contain 'm' dates with uniform infiltration values where infiltration is
    expressed in m/s. The class reads the file :

    date time1
    infiltration1
    date time2
    infiltration2
    ...
    date timem
    infiltrationm

SourceFiles
    uniformInfiltrationEventFile.C

Author
    Pierre Horgue

\*---------------------------------------------------------------------------*/

#ifndef uniformInfiltrationEventFile_H
#define uniformInfiltrationEventFile_H

#include "fvMesh.H"
#include "eventFile.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class fluiduniformInfiltrationEventFile Declaration
    \*---------------------------------------------------------------------------*/

class uniformInfiltrationEventFile
:
    public eventFile
{

public:

    // Constructors    
    //- Copy constructor
    uniformInfiltrationEventFile
    (
        const uniformInfiltrationEventFile& fileToCopy
    );

    //- Construct from Istream
    uniformInfiltrationEventFile
    (
        const word& fileName
    );
    
    // Selectors

    //- Destructor
    ~uniformInfiltrationEventFile();

    // Member Functions
    
    //- Access function
    const scalar& currentValue() {return currentValues_[0];}
    const scalar& oldValue() {return oldValues_[0];}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
