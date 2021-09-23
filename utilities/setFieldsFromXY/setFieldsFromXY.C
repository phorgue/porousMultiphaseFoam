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
    setFieldsFromXY

Description
    set a given field by performing linear interpolation using
    reference values given as table with x/y/values.
    The file should be in the form

    x1 y1 value1
    x2 y2 value2
    ...
    xN yN valueN

Usage
    setPermeabilityFieldFromXY -file inputFile_file -field fieldName
    Options are:
    -folder can be used to specify time folder
    -offset to add value to interpolated data
    -npoints to specify number of points (1: closest point, 0: number of faces
    or each cell)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "XYfile.H"

int main(int argc, char *argv[])
{

    argList::addOption("file","fileName","specify the input file");
    argList::addOption("field","fieldName","specify the output file");
    argList::addOption("folder","constant","specify the folder");
    argList::addOption("offset","0","add offset to interpolated value");
    argList::addOption("npoints","3","number of points used for interpolation (0 => number of cell faces");

    Foam::argList args(argc,argv); 

    word nameXY = "default";
    if (args.found("file"))
    {
        nameXY = args.opt("file");
    }
    else
    {
        FatalErrorIn("setFieldsFromXY.C")
            << "no input file specified, use option -file"
            << exit(FatalError);
    }

    word nameField = "default";
    if (args.found("field"))
    {
        nameField = args.opt("field");
    }

    else
    {
        FatalErrorIn("setFieldsFromXY.C")
            << "no field specified, use option -field"
            << exit(FatalError);
    }

    scalar offset = args.opt("offset",0.);
    label npoints = args.opt("npoints",3);

    #include "createTime.H"
    #include "createMesh.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- read the XY file
    XYfile sourceFile(nameXY);

    word fileDir = "constant";
    if (args.found("folder"))
    {
        fileDir = args.opt("folder");
    }

    volScalarField outputFile
        (
            IOobject
            (
                nameField,
                fileDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

    forAll(outputFile,celli)
    {
        if (npoints == 0)
        {
            outputFile[celli] = sourceFile.interpolate(mesh.C()[celli],mesh.cells()[celli].size()-2) + offset;
        }
        else
        {
            outputFile[celli] = sourceFile.interpolate(mesh.C()[celli],npoints) + offset;
        }
    }

    outputFile.write();

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
