/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    setBoundaryHeadPressure

Description
    Utility to set up Head Pressure (3D) or potential (2D) on a specified patch.
    Can be used to set up water level on lateral boundaries for (2D)groundwaterFoam

Usage
    1) for uniform head pressure (groundwaterFoam) :

      setBoundaryHeadPressure -patch patchName -value 50.3

    2) for stl dependent head pressure (groundwaterFoam)

      setBoundaryHeadPressure -patch patchName -file stl_file

    3) for stl dependent potential (groundwater2DFoam) :

      setBoundaryHeadPressure -patch patchName -file stl_file -version 2D

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurfaceMesh.H"

int main(int argc, char *argv[])
{
    argList::addOption("patch","patchName","specify the patch to set head pressure");
    argList::addOption("file","fileName","specify the STL file");
    argList::addOption("value","0","uniform potential value");
    argList::addOption("version","3D","2D (h) or 3D (potential)");

    Foam::argList args(argc,argv); 

    if (!args.optionFound("patch"))
    {
        FatalError << "no patch specified" 
            << nl << " use option -patch"
            << exit(FatalError);
    }

    #include "createTime.H"
    #include "createMesh.H"

    word version("3D");
    if (args.optionFound("version"))
    {
        version = args.option("version");
    }

    if(version == "3D")
    {
        #include "version3D.H"
    }
    else if (version == "2D")
    {
        #include "version2D.H"
    }
    else
    {
        FatalErrorIn("in setBoundaryHeadPressure.C")
            << "version not known : " << version
                << abort(FatalError);
    }

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
