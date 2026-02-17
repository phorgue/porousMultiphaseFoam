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
    Foam::flowModels::nonLinearEqn

Description
    Abtract class for non-linear equations

SourceFiles
    nonLinearEqn.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinearEqn_H
#define nonLinearEqn_H

#include "fvMesh.H"

namespace Foam
{

class nonLinearEqn
{
protected:
    const fvMesh& mesh_;

public:

    //- Constructor
    explicit nonLinearEqn(const fvMesh &mesh);

    //- Destructor
    virtual ~nonLinearEqn() = default;

    //- Access function
    const fvMesh& mesh() const { return mesh_; }

    //- Virtual functions
    virtual void updateProperties(label methodID) = 0;
    virtual Tuple2<scalar> solve(scalar tol, label methodID) = 0;
    virtual void info() = 0;

};


}

#endif

// ************************************************************************* //