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

Application
    setBoundaryHeadWaterLevel

Description
    Utility to set up head Pressure (3D) or potential (2D) on a specified patch
    Can be used to set up water level for all groundwater solvers.

Usage
    1) to impose uniform head pressure (unsaturated solvers):

        setBoundaryWaterLevel -patch patchName -value 50.3

    2) to impose uniform potential (saturated solvers):

        setBoundaryWaterLevel -patch patchName -value 50.3 -field potential

    3) initialize pressure head with DEM  file (unsaturated solvers):

        setBoundaryWaterLevel -patch patchName -DEM XYZ_file

    4) initialize potential with DEM  file (saturated solvers):

        setBoundaryWaterLevel -patch patchName -DEM XYZ_file -field potential

    5)  initialize pressure head with STL  file (unsaturated solvers):

        setBoundaryWaterLevel -patch patchName -STL stl_file

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "DEMfile.H"
#include "triSurfaceMesh.H"

int main(int argc, char *argv[])
{
    argList::addOption("patch","patchName","specify the patch to set head pressure");
    argList::addOption("STL","fileName","specify the STL file");
    argList::addOption("DEM","fileName","specify the DEM file");
    argList::addOption("value","0","uniform potential value");
    argList::addOption("field","h","h or potential");
    argList::addOption("threshold","0","minimum height for points to look in STL file");
    argList::addOption("offset","0","specify the constant offset from the STL/MNT file");
    
    Foam::argList args(argc,argv); 

    if (!args.optionFound("patch"))
    {
        FatalError << "no patch specified" 
            << nl << " use option -patch"
            << exit(FatalError);
    }

    #include "createTime.H"
    #include "createMesh.H"

    word field("h");
    if (args.optionFound("field"))
    {
        field = args.option("field");
    }

    if(field == "h")
    {
        #include "pressureheadField.H"
    }
    else if (field == "potential")
    {
        #include "potentialField.H"
    }
    else
    {
        FatalErrorIn("in setBoundaryHeadPressure.C")
            << "field not known : " << field
                << abort(FatalError);
    }

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
