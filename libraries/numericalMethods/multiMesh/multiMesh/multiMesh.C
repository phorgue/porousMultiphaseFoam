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

#include "multiMesh.H"
#include "dynamicRefineFvMesh.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMesh, 0);
    defineRunTimeSelectionTable(multiMesh, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiMesh::multiMesh (dynamicFvMesh& coarseMesh)
:
    coarseMesh_(coarseMesh)
{
    //- Protect patch cells from refinement
    if (coarseMesh_.dynamic()) {
        if (isA<dynamicRefineFvMesh>(coarseMesh_)) {
            DynamicList<label> boundary_protected_cells(0);
            const polyBoundaryMesh& patches = coarseMesh_.boundaryMesh();
            forAll(patches, patchi) {
                if (not(isA<processorPolyPatch>(patches[patchi])) &&
                    not(isA<symmetryPlanePolyPatch>(patches[patchi])) ){
                    boundary_protected_cells.append(coarseMesh_.boundary()[patchi].faceCells());
                }
            }
            refCast<dynamicRefineFvMesh>(coarseMesh_).protectedCell() = bitSet(boundary_protected_cells);
        }
    }
}

// ************************************************************************* //
