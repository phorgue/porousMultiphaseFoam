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
    groundwaterTransport2DFoam

Description
    Transient solver for free-surface flow in porous media with tracer
    transport

\*---------------------------------------------------------------------------*/

#include "PMFversion.H"
#include "fvCFD.H"
#include "fixedValueFvPatchField.H"
//#include "multiscalarMixture.H"
#include "infiltrationEventFile.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "patchEventFile.H"
#include "eventInfiltration.H"
#include "eventFlux.H"
#include "DupuitDarcyEqn.H"
#include "porousMediumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    Foam::argList args(argc, argv);
    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }
    PMFversion solverV;
    bool steady = false;

    Info << "Create time\n" << Foam::endl;
    Time runTime(Time::controlDictName, args);

    #include "createMesh.H"
    Info << "Reading transportProperties" << endl;
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
    autoPtr<incompressiblePhase> fluidPhase = incompressiblePhase::New(mesh, transportProperties, "");
    //- Porous medium model
    porousMediumModel pmModel(mesh, transportProperties, dimLength/dimTime);
    //- Dupuit-Darcy's equation
    flowModels::DupuitDarcyEqn ddEqn(mesh, transportProperties, pmModel, fluidPhase(), steady);
    //- Transport model
    Info << "Reading composition" << endl;
    autoPtr<porousMediumTransportModel> pmTransportModel =
            porousMediumTransportModel::New("", mesh, transportProperties);
    multiscalarMixture& composition = pmTransportModel->composition();

    //- create water-source/infiltration/tracer-source events
    autoPtr<infiltrationEventFile> infiltrationEvent = infiltrationEventFile::New("infiltrationEventFile", transportProperties);
    infiltrationEvent->init(runTime, ddEqn.potential().name(), mesh, ddEqn.infiltration());
    autoPtr<sourceEventFile> waterSourceEvent = sourceEventFile::New("sourceEventFileWater", transportProperties);
    waterSourceEvent->init(runTime, ddEqn.potential().name(), mesh, dimLength/dimTime);
    List<sourceEventFile*>& tracerSourceEventList = pmTransportModel->sourceEventList();
    forAll(tracerSourceEventList,sourceEventi) tracerSourceEventList[sourceEventi]->init(runTime);

    //- create time manager
    List<sourceEventFile*> fullSourceEventList;
    fullSourceEventList.append(waterSourceEvent.get());
    fullSourceEventList.append(tracerSourceEventList);
    List<infiltrationEventFile*> infiltrationEventList;
    infiltrationEventList.append(infiltrationEvent.get());
    multiDtManager MDTM(runTime, fullSourceEventList, infiltrationEventList);
    MDTM.addField(ddEqn.hwater(), ddEqn.dryCellIDList());
    forAll(composition.Y(), speciesi) MDTM.addField(composition.Y()[speciesi], ddEqn.dryCellIDList());

    autoPtr<outputEventFile> outputEvent = outputEventFile::New(runTime, mesh, ddEqn.zScale());
    outputEvent->addField(ddEqn.hwater(), ddEqn.phi(), pmModel.eps(), "waterMassBalance.csv");
    outputEvent->addSourceTerm("fixedPoints", ddEqn.flowInOutFixedPoints());
    outputEvent->addSourceTerm("seepage", ddEqn.flowOutSeepage());
    outputEvent->addField(ddEqn.potential(), ddEqn.phi());
    forAll(composition.Y(), speciei) {
        outputEvent->addField(composition.Y()[speciei], ddEqn.phihwater(), pmModel.eps(), ddEqn.hwater(), composition.R(speciei), composition.Y()[speciei].name()+"massBalance.csv");
    //    outputEvent->addSourceTerm("seepage", ddEqn.outflowSeepageTracer[speciei]);
    }
    outputEvent->init();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        if (infiltrationEvent->isPresent()) infiltrationEvent->updateIndex(runTime.timeOutputValue());
        if (waterSourceEvent->isPresent()) waterSourceEvent->updateIndex(runTime.timeOutputValue());
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
        forAll(tracerSourceEventList,sourcei) tracerSourceEventList[sourcei]->updateIndex(runTime.timeOutputValue());

        MDTM.updateDt();

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        //- Update water infiltration + source
        if (infiltrationEvent->isPresent()) infiltrationEvent->updateInfiltration(runTime, ddEqn.infiltration());
        if (waterSourceEvent->isPresent()) {
            waterSourceEvent->updateValue(runTime);
            pmModel.sourceTerm() = waterSourceEvent->dtValuesAsField();
        }

        //- Solve Dupuit-Darcy potential equation
        ddEqn.solve();

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

//    if (cumulativeWaterAdded > 0) Info << "Cumulated water added = " << cumulativeWaterAdded << " m3, equivalent height = " << cumulativeWaterAdded*zScale/gSum(mesh.V()) << " m" << nl << endl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
