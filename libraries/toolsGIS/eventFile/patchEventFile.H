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
    Foam::patchEventFile

Description
    PatchEvent file handler/reader for the porousMultiphaseFoam toolbox which
    contain 'm' dates with 'n' infiltration points where infiltration is
    expressed in m3/s. The class reads the file :

    date time1
    patchName1 infiltration11
    patchName1 infiltration12
    ...
    patchNameN infiltration1N
    date time2
    patchName1 infiltration21
    patchName2 infiltration22
    ...
    date time
    ...
    patchNameN infiltrationmN

SourceFiles
    patchEventFile.C

Author
    Pierre Horgue

\*---------------------------------------------------------------------------*/

#ifndef patchEventFile_H
#define patchEventFile_H

#include "fvMesh.H"
#include "eventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class fluidpatchEventFile Declaration
    \*---------------------------------------------------------------------------*/

class patchEventFile
:
    public eventFile
{

protected:

    label npatches_;
    wordList patchNameList_;

public:

    // Constructors    
    //- Copy constructor
    patchEventFile
    (
        const patchEventFile& fileToCopy
    );

    //- Construct from Istream
    patchEventFile
    (
        const word& fileName,
        bool display = true
    );
    
    // Selectors

    //- Destructor
    ~patchEventFile();

    // Member Functions
    
    //- Access function
    const label& npatches() const {return npatches_;}
    const wordList& patchNameList() const {return patchNameList_;}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
