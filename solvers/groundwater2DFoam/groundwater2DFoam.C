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

Application
    groundwater2DFoam

Description
    Transient solver for free-surface flow in porous media

Developer
    P. Horgue

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "harmonic.H"
#include "fixedValueFvPatchField.H"
#include "MNTfile.H"
#include "infiltrationEventFile.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "EulerD3dt3Scheme.H"
#include "EulerD2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "readFixedPoints.H"
    #include "readEvent.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        if (infiltrationEventIsPresent) infiltrationEvent.updateIndex(runTime.timeOutputValue());
        if (waterSourceEventIsPresent) waterSourceEvent.updateIndex(runTime.timeOutputValue());
        #include "setDeltaT.H"

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        //- Update infiltration term
        #include "computeInfiltration.H"

        //- Solve potential equation
        #include "potentialEqn.H"

        //- Water mass balance computation
        #include "waterMassBalance.H"

        #include "eventWrite.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    if (cumulativeWaterAdded > 0) Info << "Cumulated water added = " << cumulativeWaterAdded << " m3, equivalent height = " << cumulativeWaterAdded*zScale/sum(mesh.V()).value() << " m" << nl << endl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
