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
    porouScalarTransport2DFoam

Description
    Solves the transport equation for a passive scalar in porous media 
    with dispersion coefficient model considering a 2D groundwater model 
    Velocity field should be pre-computed using 2D solver, groundater2DFoam
    or steadyGroundwater2DFoam

\*---------------------------------------------------------------------------*/

#include "PMFversion.H"
#include "fvCFD.H"
#include "outputEventFile.H"
#include "eventFlux.H"
#include "multiDtManager.H"
#include "frozenDupuitDarcyEqn.H"
#include "porousMediumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList args(argc, argv);
    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }
    PMFversion solverV;

    Info<< "Create time\n" << Foam::endl;
    Time runTime(Time::controlDictName, args);

    #include "createMesh.H"
    IOdictionary transportProperties
    (
    IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<patchEventFile*> patchEventList;
    eventFlux::setEventFileRegistry(&patchEventList, "C");

    //- fluid phase model
    autoPtr<fluidPhase> fluidPhase = fluidPhase::New(mesh, transportProperties, "");
    //- Porous medium model
    porousMediumModel pmModel(mesh, transportProperties, dimLength/dimTime);
    //- Frozen flow-field (Dupuit-Darcy's equation)
    flowModels::frozenDupuitDarcyEqn ddEqn(mesh, transportProperties, pmModel, fluidPhase.ref());
    //- Transport model
    Info << "Reading composition" << endl;
    autoPtr<porousMediumTransportModel> pmTransportModel =
            porousMediumTransportModel::New("", mesh, transportProperties);
    multiscalarMixture& composition = pmTransportModel->composition();

    //- create tracer-source events
    List<sourceEventFile*>& tracerSourceEventList = pmTransportModel->sourceEventList();
    forAll(tracerSourceEventList,sourceEventi) tracerSourceEventList[sourceEventi]->init(runTime);
    forAll(patchEventList,patchEventi) patchEventList[patchEventi]->init(runTime);

    //- Create timestep manager
    multiDtManager MDTM(runTime, tracerSourceEventList, patchEventList);
    forAll(composition.Y(), speciesi) MDTM.addField(composition.Y()[speciesi]);

    //- Output events
    autoPtr<outputEventFile> outputEvent = outputEventFile::New(runTime, mesh, ddEqn.zScale());
    forAll(composition.Y(), speciei) {
        outputEvent->addField(composition.Y()[speciei], ddEqn.phihwater(), pmModel.eps(), ddEqn.hwater(), composition.R(speciei), composition.Y()[speciei].name()+"massBalance.csv");
     //   outputEvent->addSourceTerm("seepage", outflowSeepageTracer[speciei]);
    }
    outputEvent->init();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.run())
    {
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
        forAll(tracerSourceEventList,sourceEventi) tracerSourceEventList[sourceEventi]->updateIndex(runTime.timeOutputValue());

        MDTM.updateDt();

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        //- Solve transport equation
        forAll(tracerSourceEventList, sourcei) tracerSourceEventList[sourcei]->updateValue(runTime);
        pmTransportModel->solveTransport(fluidPhase->U(),
                                         ddEqn.phihwater(),
                                         pmModel.eps(),
                                         ddEqn.hwater(),
                                         ddEqn.seepage(),
                                         ddEqn.zScale());

        MDTM.updateAllDerivatives();

        //- Write
        outputEvent->write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
