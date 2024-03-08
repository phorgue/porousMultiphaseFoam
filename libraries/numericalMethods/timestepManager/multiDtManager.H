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
    Foam::fv::multiDtManager

Description
    Class for manage timestep using multiple timestepManagerTruncation (for each field)

SourceFiles
    multiDtManager.C

\*---------------------------------------------------------------------------*/

#ifndef multiDtManager_H
#define multiDtManager_H

#include "timestepManagerTruncation.H"
#include "timestepManagerIterative.H"
#include "sourceEventFile.H"
#include "patchEventFile.H"
#include "infiltrationEventFile.H"

namespace Foam
{

class multiDtManager
{

private:
    Time& runTime_;

    bool adjustTimeStep_;
    scalar maxDeltaT_;
    bool eventTimeTracking_;

    PtrList<timestepManagerTruncation> dtManagerT_;
    PtrList<timestepManagerIterative> dtManagerI_;
    const List<sourceEventFile*> sourceEventList_;
    const List<patchEventFile*> patchEventList_;
    const List<infiltrationEventFile*> infiltrationEventList_;

public:
    // Constructors

    //- Construct from field
    multiDtManager(Time&, const List<sourceEventFile*>&, const List<patchEventFile*>&);
    multiDtManager(Time&, const List<sourceEventFile*>&, const List<infiltrationEventFile*>&);

    //- Destructor
    virtual ~multiDtManager();

    //- Access
    timestepManagerTruncation& dtManagerT(label index){return dtManagerT_[index];};
    timestepManagerIterative& dtManagerI(label index){return dtManagerI_[index];};
    const bool& adjustTimeStep() {return adjustTimeStep_;};
    
    //- Functions
    void addField(const volScalarField&, const labelList* dryCells = nullptr);
    void addField(const volScalarField&, const labelList& dryCells);
    void addIterativeAlgorithm(const volScalarField&, const word&);
    void updateDt();
    void updateAllDerivatives();
    void adjustDeltaTUsingEvent();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
