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

#include "eventFile.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eventFile::eventFile
(
    const eventFile& eventFileToCopy
)
    :  
    name_(eventFileToCopy.name()),
    ndates_(eventFileToCopy.ndates()),
    dates_(eventFileToCopy.dates()),
    datas_(eventFileToCopy.datas_),
    currentValues_(eventFileToCopy.currentValues_),
    oldValues_(eventFileToCopy.oldValues_),
    iterator_(eventFileToCopy.iterator_)
{
}

Foam::eventFile::eventFile
(
    const word& fileName
)
    :
    name_(fileName),
    ndates_(0),
    dates_(),
    datas_(),
    currentValues_(),
    oldValues_(),
    iterator_(-1)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eventFile::~eventFile()
{}

// * * * * * * * * * * * * * * * * Members  * * * * * * * * * * * * * * * //

Foam::scalar Foam::eventFile::currentEventStartTime() const
{
    if (iterator_ == -1)
    {
        return 0;
    }
    else
    {
        return dates_[iterator_];
    }
}

const Foam::scalar& Foam::eventFile::currentEventEndTime() const
{
    if (iterator_ < ndates_-1)
    {
        return dates_[iterator_+1];
    }
    else
    {
        return GREAT;
    }
}

void Foam::eventFile::update(const scalar& currentTime)
{
    storeOldValues();
    if (currentEventEndTime() <= currentTime)
    {
        iterator_++;
        while ((currentEventEndTime() <= currentTime) && (iterator_ < ndates_-2))
        {
            iterator_++;
        }
    }
    if (currentTime < dates_[0])
    {
        currentValues_ = 0.0;
    }
    else if (iterator_ < ndates_-1)
    {
        scalar interpolateFactor_ = (currentTime - dates_[iterator_]) / (dates_[iterator_+1] - dates_[iterator_]);
        forAll(currentValues_,id)
        {
            currentValues_[id] = (1.0 - interpolateFactor_) * datas_[iterator_][id] + interpolateFactor_ * datas_[iterator_+1][id];
        }
    }
    else
    {
        currentValues_ = 0.0;
    }
}

void Foam::eventFile::addIntermediateTimeSteps(const scalar& smallDeltaT)
{
    RectangularMatrix<scalar> oldDatas = datas_;
    scalarList oldDates = dates_;
    datas_.setSize((ndates_-2)*3+2,datas_.n());
    dates_.setSize((ndates_-2)*3+2);
    dates_[0] = oldDates[0];
    for (label datei=1;datei<ndates_-1;datei++)
    {
        dates_[datei*3-2] = oldDates[datei]-smallDeltaT;
        dates_[datei*3-1] = oldDates[datei];
        dates_[datei*3] = oldDates[datei]+smallDeltaT;
    }
    dates_[(ndates_-2)*3+1] = oldDates[ndates_-1];
    for(label columni=0;columni<datas_.n();columni++)
    {
        datas_[0][columni] = oldDatas[0][columni];
        for(label datei=1;datei<ndates_-1;datei++)
        {
            datas_[datei*3-2][columni] = (oldDatas[datei-1][columni]+oldDatas[datei][columni])/2;
            datas_[datei*3-1][columni] = (oldDatas[datei-1][columni]+oldDatas[datei][columni])/2;
            datas_[datei*3][columni] = oldDatas[datei][columni];
        }
        datas_[(ndates_-2)*3+1][columni] = oldDatas[ndates_-1][columni];
    }
    ndates_ = (ndates_-2)*3+2;
}
