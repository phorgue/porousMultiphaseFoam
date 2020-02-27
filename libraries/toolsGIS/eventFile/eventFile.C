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
#include "EulerDdtScheme.H"
#include "steadyStateDdtScheme.H"
#include "backwardDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"

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
    oldOldValues_(eventFileToCopy.oldOldValues_),
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
    oldOldValues_(),
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

void Foam::eventFile::updateIndex(const scalar& currentTime)
{
    if  (iterator_ < ndates_-1)
    {
        while (currentEventEndTime() <= currentTime)
        {
            iterator_++;
            if (iterator_ == ndates_-1) break;
        }
    }
}

void Foam::eventFile::updateValue(const TimeState& runTime)
{
    storeOldValues();
    if (runTime.timeOutputValue() < dates_[0])
    {
        currentValues_ = 0.0;
    }
    else if (iterator_ == -1)
    {
        scalar dt2 = runTime.timeOutputValue() - dates_[0];
        forAll(currentValues_,id)
        {
            scalar value2 = datas_[0][id] + dt2 * (datas_[1][id]-datas_[0][id])/(dates_[1]-dates_[0]);
            currentValues_[id] = dt2 * value2 / runTime.deltaTValue();
        }
    }
    else if (iterator_ < ndates_-1)
    {
        if (runTime.timeOutputValue() <= dates_[iterator_+1])
        {
            //- T and T+deltaT in the same event
            scalar interpolateFactor = (runTime.timeOutputValue() - runTime.deltaTValue()/2. - dates_[iterator_]) / (dates_[iterator_+1] - dates_[iterator_]);
            forAll(currentValues_,id)
            {
                currentValues_[id] = (1.0 - interpolateFactor) * datas_[iterator_][id] + interpolateFactor * datas_[iterator_+1][id];
            }
        }
        else
        {
            //- T and T+deltaT in different events
            scalar dt1 = dates_[iterator_+1] - (runTime.timeOutputValue()-runTime.deltaTValue());
            scalarList value1(currentValues_.size(),0.);
            forAll(currentValues_,id) value1[id] = datas_[iterator_+1][id] + dt1 * (datas_[iterator_][id]-datas_[iterator_+1][id])/(dates_[iterator_]-dates_[iterator_+1]);
            scalar dt2 = 0;
            scalarList value2(currentValues_.size(),0.);
            if (iterator_ < ndates_-2)
            {
                //- handling same dates with different values (heavy-side functions)
                label iteratorNext = iterator_+1;
                if (dates_[iteratorNext] == dates_[iteratorNext+1]) iteratorNext++;

                if (iteratorNext == ndates_-1) FatalErrorIn("eventFile.C") << "event file : " << this->name() << " finished by two same dates, remove the last one" << abort(FatalError);

                scalar dt2 = runTime.timeOutputValue() - dates_[iteratorNext];
                forAll(currentValues_,id) value2[id] = datas_[iteratorNext][id] + dt2 * (datas_[iteratorNext+1][id]-datas_[iteratorNext][id])/(dates_[iteratorNext+1]-dates_[iteratorNext]);
            }
            forAll(currentValues_,id)
            {
                currentValues_[id] = (dt1 * value1[id] + dt2 * value2[id]) / runTime.deltaTValue();
            }
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


void Foam::eventFile::setTimeScheme(const word& dtFieldName, const fvMesh& mesh)
{
    ddtScheme_ = fv::ddtScheme<scalar>::New
    ( 
        mesh,
        mesh.ddtScheme("ddt(" + dtFieldName + ')')
    );

    mesh_ = &mesh;
    onMeshChanged();
}

Foam::scalar Foam::eventFile::dtValue(const label& id) const
{
    if(ddtScheme_.empty())
    {
        FatalErrorIn("eventFile.C")
            << "You must call setTimeScheme(...) before being able to use dtValue(s)()"
            << abort(FatalError);
    }

    const fv::ddtScheme<scalar>& scheme = ddtScheme_.ref();

    using Euler = fv::EulerDdtScheme<scalar>;
    using steadyState = fv::steadyStateDdtScheme<scalar>;
    using backward = fv::backwardDdtScheme<scalar>;
    using CrankNicolson = fv::CrankNicolsonDdtScheme<scalar>;

    if (dynamic_cast<const backward*>(&scheme))
    {
        scalar deltaT = mesh_->time().deltaT().value();
        scalar deltaT0 = mesh_->time().deltaT0().value();
        
        scalar coefft0_00 = deltaT/(deltaT + deltaT0);
        scalar coefftn_0 = 1 + coefft0_00;

        return coefftn_0*this->currentValue(id) - coefft0_00*this->oldValue(id);
    }
    else if (const auto CNscheme = dynamic_cast<const CrankNicolson*>(&scheme))
    {
        const auto& ocCoeff = CNscheme->ocCoeff();
        return (1 + ocCoeff)*this->currentValue(id) - ocCoeff*(1+ocCoeff) * this->oldValue(id) + ocCoeff*ocCoeff*this->oldOldValue(id);
    }
    else if (dynamic_cast<const Euler*>(&scheme) || dynamic_cast<const steadyState*>(&scheme))
    {
        return this->currentValue(id);
    }

    FatalErrorIn("eventFile.C")
        << "ddtScheme " << scheme.type() << " unsupported"
        << abort(FatalError);
    return 0;
}

Foam::scalarList Foam::eventFile::dtValues() const
{

    scalarList ret(currentValues().size());
    
    forAll(ret, id)
    {
        ret[id] = dtValue(id);
    }

    return ret;
}

