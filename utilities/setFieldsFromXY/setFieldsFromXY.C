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
    setPermeabilityFieldFromXY

Description
    set the permeability field by performing linear interpolation using 
    permeability reference values (as table with x,y,K)
    For N points, the file should be in the form
    N 3
    (
    (x1 y1 K1)
    (x2 y2 K2)
    ...
    (xN yN kN)
    )

Usage
    setPermeabilityFieldFromXY -file inputFile_file

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    argList::addOption("fileIn","fileName","specify the input file");
    argList::addOption("fileOut","fileName","specify the output file");   
    argList::addOption("folder","constant","specify the folder");

    Foam::argList args(argc,argv); 

    if (!args.optionFound("fileIn"))
    {
        FatalError << "no input file specified" 
            << nl << " use option -fileIn"
            << exit(FatalError);
    }
    else if (!args.optionFound("fileOut"))
    {
        FatalError << "no output file specified" 
            << nl << " use option -fileOut"
            << exit(FatalError);
    }
    
    #include "createTime.H"
    #include "createMesh.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    //- read inputFile permeability matrix
    RectangularMatrix<scalar> inputFile(IFstream(args.option("fileIn"))());

    word fileDir = "constant";
    if (args.optionFound("folder"))
    {
        fileDir = args.option("folder"); 
    
    }

    volScalarField outputFile
        (
            IOobject
            (
                args.option("fileOut"),
                fileDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

    forAll(outputFile,celli)
    {   
        label id1=-1;
        label id2=-1;
        label id3=-1;     
        scalar dist1 = VGREAT;
        scalar dist2 = VGREAT;
        scalar dist3 = VGREAT;

        for(label i=0;i<inputFile.m();i++)
        {
            scalar dist = Foam::sqrt(pow(inputFile[i][0]-mesh.C()[celli].x(),2)+pow(inputFile[i][1]-mesh.C()[celli].y(),2));
            if ( dist < dist1)
            {
                //Info <<  "1 *** " << dist << " " << dist1 <<  " " << dist2 <<  " " << dist3 << endl;
                id3 = id2;
                dist3 = dist2;
                id2 = id1;
                dist2 = dist1;
                id1 = i;
                dist1 = dist;
            }
            else if ( dist < dist2)
            {
                //         Info <<  "2 *** " << dist << " " << dist1 <<  " " << dist2 <<  " " << dist3 << endl;      
                id3 = id2;
                dist3 = dist2;
                id2 = i;
                dist2 = dist;                
            }
            else if ( dist < dist3)
            {
                //                Info <<  "3 *** " << dist << " " << dist1 <<  " " << dist2 <<  " " << dist3 << endl;
                id3 = i;
                dist3 = dist;                
            }
        }

        
        if ( (id1 == -1) || (id2 == -1) || (id3 == -1))
        {
            Info << nl << "Error : three point are not found for interpolation"
                << nl << id1 << " / " << id2 << " / " << id3 << endl;
        }
       
        outputFile[celli] = ( dist1*inputFile[id1][2] + dist2*inputFile[id2][2] + dist3*inputFile[id3][2] ) / (dist1+dist2+dist3);

    }

    outputFile.write();
    
    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
