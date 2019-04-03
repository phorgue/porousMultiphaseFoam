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

protected:

    word name_;
    label ndates_;
    scalarList dates_;
    RectangularMatrix<scalar> datas_;
    scalarList currentValues_;
    scalarList oldValues_;
    label iterator_;

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
    const scalarList& dates() const {return dates_;}
    const RectangularMatrix<scalar>& datas() const {return datas_;}
    const scalar& lastValue(const label& id) const {return datas_[iterator_][id];}
    const scalar& oldValue(const label& id) const {return oldValues_[id];}
    const scalar& currentValue(const label& id) const {return currentValues_[id];}
    const scalar& currentEventStartTime() const;
    const scalar& currentEventEndTime() const;
    void storeOldValues() {oldValues_ = currentValues_;}
    void update(const scalar& currentTime);
    void addIntermediateTimeSteps(const scalar& smallDT);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
