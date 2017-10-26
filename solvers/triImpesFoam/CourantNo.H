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
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNuma = 0.0;
scalar meanCoNumw = 0.0;
scalar meanCoNumo = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhia
    (
        fvc::surfaceSum(mag(phia))().internalField()
    );
    scalarField sumPhiw
    (
        fvc::surfaceSum(mag(phiw))().internalField()
    );
    scalarField sumPhio
    (
        fvc::surfaceSum(mag(phio))().internalField()
    );
    
    CoNum = max(
        0.5*gMax(sumPhia/mesh.V().field())*runTime.deltaTValue(),
        max(
            0.5*gMax(sumPhiw/mesh.V().field())*runTime.deltaTValue(),
            0.5*gMax(sumPhio/mesh.V().field())*runTime.deltaTValue()
        )
    );

    meanCoNuma = 
        0.5*(gSum(sumPhia)/gSum(mesh.V().field()))*runTime.deltaTValue();
    meanCoNumw = 
        0.5*(gSum(sumPhiw)/gSum(mesh.V().field()))*runTime.deltaTValue();
    meanCoNumo = 
        0.5*(gSum(sumPhio)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNuma << "," << meanCoNumw << "," << meanCoNumo << " max: " << CoNum << endl;

maxDeltaTFact = maxCo/(CoNum + SMALL);

// ************************************************************************* //
