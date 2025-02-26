/*---------------------------------------------------------------------------* \
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
    groundwater2DFoam

Description
    Transient solver for free-surface flow in porous media

\*---------------------------------------------------------------------------*/

#include "PMFversion.H"
#include "fvCFD.H"
#include "fixedValueFvPatchField.H"
#include "infiltrationEventFile.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "multiDtManager.H"
#include "DupuitDarcyEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addBoolOption("steady", "to run steady flow simulation");
    Foam::argList args(argc, argv);

    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }
    PMFversion solverV;
    bool steady = args.found("steady");

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

    //- fluid phase model
    autoPtr<incompressiblePhase> fluidPhase = incompressiblePhase::New(mesh, transportProperties, "");
    //- 2D porous medium model
    porousMediumModel pmModel(mesh, transportProperties, dimLength/dimTime);
    //- Richards' equation
    flowModels::DupuitDarcyEqn ddEqn(mesh, transportProperties, pmModel, fluidPhase(), steady);

    const dictionary& residualControl = mesh.solutionDict().subOrEmptyDict("residualControl");
    const auto residualPotential = residualControl.lookupOrDefault<scalar>("potential", 0);
    if (steady && residualPotential ==0)
    {
        FatalErrorIn("readTimeControls.h") << "residualControl.potential should be specified in system/fvSolution" << abort(FatalError);
    }

    //- create source/infiltration events
    autoPtr<infiltrationEventFile> infiltrationEvent = infiltrationEventFile::New("infiltrationEventFile", transportProperties);
    infiltrationEvent->init(runTime, ddEqn.potential().name(), mesh, ddEqn.infiltration());
    autoPtr<sourceEventFile> sourceEvent = sourceEventFile::New("sourceEventFileWater", transportProperties);
    sourceEvent->init(runTime, ddEqn.potential().name(), mesh, dimLength/dimTime);

    //- create time manager
    List<sourceEventFile*> sourceEventList;
    sourceEventList.append(sourceEvent.get());
    List<infiltrationEventFile*> infiltrationEventList;
    infiltrationEventList.append(infiltrationEvent.get());
    multiDtManager MDTM(runTime, sourceEventList, infiltrationEventList);
    MDTM.addField(ddEqn.potential(), ddEqn.dryCellIDList());

    autoPtr<outputEventFile> outputEvent = outputEventFile::New(runTime, mesh, ddEqn.zScale());
    outputEvent->addField(ddEqn.hwater(), ddEqn.phi(), pmModel.eps(), "waterMassBalance.csv");
    outputEvent->addSourceTerm("fixedPoints", ddEqn.flowInOutFixedPoints());
    outputEvent->addSourceTerm("seepage", ddEqn.flowOutSeepage());
    outputEvent->addField(ddEqn.potential(), ddEqn.phi());
    outputEvent->init();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        if (!steady)
        {
            if (infiltrationEvent->isPresent()) infiltrationEvent->updateIndex(runTime.timeOutputValue());
            if (sourceEvent->isPresent()) sourceEvent->updateIndex(runTime.timeOutputValue());
            MDTM.updateDt();
        }

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        //- Update infiltration term
        if (!steady)
        {
            if (infiltrationEvent->isPresent()) infiltrationEvent->updateInfiltration(runTime, ddEqn.infiltration());

           if (sourceEvent->isPresent()) {
                sourceEvent->updateValue(runTime);
                pmModel.sourceTerm() = sourceEvent->dtValuesAsField();
            }
        }

        //- Solve Dupuit-Darcy equation
        scalar residual = ddEqn.solve();

        //- Residual computation
        if (steady)
        {
            if (residual < residualPotential) runTime.writeAndEnd();
            else runTime.write();
        }
        else
        {
            MDTM.updateAllDerivatives();
            outputEvent->write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
