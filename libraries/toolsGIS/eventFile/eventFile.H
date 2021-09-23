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
    Foam::eventFile

Description
    The is an general class used to handle/read the various event files for
    the porousMultiphaseFoam toolbox. This can be used to impose volumic or 
    mass flow rate source inside the domain or on specified patch. It can be
    used to impose time variable uniform infiltration or specify user-defined
    time outputs.

SourceFiles
    eventFile.C

\*---------------------------------------------------------------------------*/

#ifndef eventFile_H
#define eventFile_H

#include "fvMesh.H"
#include "TimeState.H"
#include "ddtScheme.H"

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
    scalarList oldOldValues_;
    label iterator_;

    tmp<fv::ddtScheme<scalar>> ddtScheme_;
    const fvMesh* mesh_ = nullptr;

    virtual void onMeshChanged() {};

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
    scalar iterator() const {return iterator_;}
    const RectangularMatrix<scalar>& datas() const {return datas_;}
    const scalar& lastValue(const label& id) const {return datas_[iterator_][id];}
    const scalar& oldValue(const label& id) const {return oldValues_[id];}
    const scalar& oldOldValue(const label& id) const {return oldOldValues_[id];}
    const scalarList& oldValues() const {return oldValues_;}
    const scalar& currentValue(const label& id) const {return currentValues_[id];}
    const scalarList& currentValues() const {return currentValues_;}
    scalar currentEventStartTime() const;
    const scalar& currentEventEndTime() const;
    void storeOldValues() {
        oldOldValues_ = oldValues_;
        oldValues_ = currentValues_;
    }
    void updateIndex(const scalar& currentTime);
    void updateValue(const TimeState& currentTime);
    void addIntermediateTimeSteps(const scalar& smallDT);

    //- Compute values from current and old considering the time discretization

    void setTimeScheme(const word& dtFieldName, const fvMesh& mesh);

    scalar dtValue(const label& id) const;
    scalarList dtValues() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
