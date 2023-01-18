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

\*---------------------------------------------------------------------------*/

#include "twophasePorousMediumModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::twophasePorousMediumModel> Foam::twophasePorousMediumModel::New
(
    const word Sname,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const autoPtr<incompressiblePhase>& phase,
    const word porousRegion
)
{
    const word modelType(transportProperties.lookupOrDefault<word>("porousMediumModel", "simplePorosity"));

    Info<< "Selecting twophasePorousMedium model => " << modelType << "\n" << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
            (
                "twophasePorousMediumModel::New(...)"
            )   << "Unknown twophasePorousMediumModel type "
                << modelType << nl << nl
                << "Valid twophasePorousMediumModels are : " << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
    }

    return autoPtr<twophasePorousMediumModel>
        (cstrIter()(Sname, mesh, transportProperties, phase, porousRegion));
}


// ************************************************************************* //
