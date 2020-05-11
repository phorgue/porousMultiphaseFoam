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

Global
    readPicardControls

Description
    Read control for Picard's iterations

\*---------------------------------------------------------------------------*/

const dictionary& PicardDict = mesh.solutionDict().subOrEmptyDict("Picard");
const scalar tolerancePicard = PicardDict.lookupOrDefault<scalar>("tolerance", 100);
const label maxIterPicard = PicardDict.lookupOrDefault<label>("maxIter",10);

Info << nl << "Picard loop control" << nl << "{"
<< nl << "    tolerance = " << tolerancePicard
<< "    maximum number of iteration = " << maxIterPicard
<< nl << "}" << endl;

if (!PicardDict.found("tolerance"))
{
    FatalErrorIn("readPicardControls.H")
        << "Missing keyword tolerance in system/fvSolution.Picard"
            << abort(FatalError);
}

//- to switch between conservative mixed formulation (by default) and head pressure (non conservative)
bool massConservative = transportProperties.lookupOrDefault<bool>("massConservative",true);
scalar massConservativeTerms = 0;
if (massConservative) massConservativeTerms = 1;
