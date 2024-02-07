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
    porousScalarTransportFoam

Description
    Solves the transport equation for a passive scalar
    in porous media with dispersion coefficient model

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidPhase.H"
#include "porousMediumTransportModel.H"
#include "multiscalarMixture.H"
#include "sourceEventFile.H"
#include "patchEventFile.H"
#include "outputEventFile.H"
#include "eventFlux.H"
#include "multiDtManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList args(argc, argv);
    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }
    #include "../headerPMF.H"

    Info<< "Create time\n" << Foam::endl;
    Time runTime(Time::controlDictName, args);

    #include "createMesh.H"
    #include "createFields.H"
    multiDtManager MDTM(runTime, tracerSourceEventList, patchEventList);
    forAll(composition.Y(), speciesi) MDTM.addField(composition.Y()[speciesi]);
    #include "readEvent.H"
    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.run())
    {
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
        forAll(tracerSourceEventList,sourceEventi) tracerSourceEventList[sourceEventi]->updateIndex(runTime.timeOutputValue());
        MDTM.updateDt();

        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "CEqn.H"
        #include "CmassBalance.H"

        #include "eventWrite.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
