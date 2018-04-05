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
    setHeadPressureFoam

Description
    Utility to set up Head Pressure on a specified patch. Can be used to
    set up water level on lateral boundaries for groundwaterFoam

    \*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurfaceMesh.H"

int main(int argc, char *argv[])
{
    argList::addOption("patch","patchName","specify the patch to set head pressure");
    argList::addOption("file","fileName","specify the STL file");
    argList::addOption("value","0","uniform potential value");
  
    Foam::argList args(argc,argv); 

    if (!args.optionFound("patch"))
    {
        FatalError << "no patch specified" 
            << nl << " use option -patch"
            << exit(FatalError);
    }

    #include "createTime.H"
    #include "createMesh.H"

    volScalarField h
        (
            IOobject
            (
                "h",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

  
    //-- Reading patch information
    word patchName = args.option("patch");
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    fvPatchScalarField& hPatch = h.boundaryFieldRef()[patchID];
    const vectorField& faces = mesh.boundary()[patchID].patch().faceCentres();
  
    //-- Compute and set head pressure
    if (args.optionFound("file"))
    {
        //- reading STL informations
        word STLfile = args.option("file");
        triSurfaceMesh potentialSTL(IOobject(STLfile,mesh));
        pointField pPoints = potentialSTL.points();
        Info << nl << "Potential fixed using STL = " << STLfile << endl;
      
        //- computing local potential
        forAll(hPatch,facei)
        {
            scalar xy_distance = GREAT;
            label id_point = -1;
            forAll(pPoints,pointi)
            {
                scalar tmp_dist = Foam::sqrt(pow(pPoints[pointi].x()-faces[facei].x(),2)+pow(pPoints[pointi].y()-faces[facei].y(),2));
                if (tmp_dist < xy_distance)
                {
                    xy_distance = tmp_dist;
                    id_point = pointi;
                }
            }
            hPatch[facei] = pPoints[id_point].z() - faces[facei].z();
        }

    }  
    else
    {
        scalar potential = args.optionLookupOrDefault<scalar>("value",0.);
        Info << nl << "Uniform potential fixed = " << potential << " m " << endl;
        forAll(faces,facei)
        {
            hPatch[facei] = (potential - faces[facei].z());
        }
    }

    h.write();

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
