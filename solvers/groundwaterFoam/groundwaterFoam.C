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
    groundwaterFoam

Description
    Transient or steady solver for Richards equation.
    A Picard loop is used for linearization.
    Permeability is isotropic (K == volScalarField)

\*---------------------------------------------------------------------------*/

#include "PMFversion.H"
#include "fvCFD.H"
#include "incompressiblePhase.H"
#include "twophasePorousMediumModel.H"
#include "fixedValueFvPatchField.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "multiDtManager.H"
#include "RichardsEqn.H"
#include "eventInfiltration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addBoolOption("steady", "to run steady flow simulation");

    Foam::argList args(argc, argv);
    bool steady = args.found("steady");
    PMFversion solverV;

    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }

    Info << "Create time\n" << Foam::endl;
    Time runTime(Time::controlDictName, args);
    bool writeResiduals(runTime.controlDict().getOrDefault<bool>("writeResiduals", false));

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
    eventInfiltration::setEventFileRegistry(&patchEventList, "h");

    //- fluid phase model
    autoPtr<incompressiblePhase> fluidPhase = incompressiblePhase::New(mesh, transportProperties, "theta");
    //- two-phase flow model
    autoPtr<twophasePorousMediumModel> pmModel =
        twophasePorousMediumModel::New("theta", mesh, transportProperties, fluidPhase);
    //- Richards' equation
    flowModels::RichardsEqn hEqn(mesh, transportProperties, pmModel(), fluidPhase(), steady);

    //- create source event for water
    autoPtr<sourceEventFile> sourceEvent = sourceEventFile::New("sourceEventFileWater", transportProperties);
    sourceEvent->init(runTime, hEqn.h().name(), mesh, pmModel->sourceTerm().dimensions());

    //- create time managers
    List<sourceEventFile*> sourceEventList;
    sourceEventList.append(sourceEvent.get());

    multiDtManager MDTM(runTime, sourceEventList, patchEventList);
    MDTM.addIterativeAlgorithm(hEqn.theta(), "Picard", steady);
    MDTM.addIterativeAlgorithm(hEqn.theta(), "Newton", steady);
    const labelList* fixedPotentialIDListPtr  = &hEqn.seepageIDList();
    MDTM.addField(hEqn.h(), fixedPotentialIDListPtr);

    //- output event
    autoPtr<outputEventFile> outputEvent = outputEventFile::New(runTime, mesh);
    outputEvent->addField(hEqn.h(), hEqn.phi());
    outputEvent->addField(hEqn.theta(), hEqn.phi(), "waterMassBalance.csv", true);
    outputEvent->init();

    OFstream residualFile("residuals.csv");

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    timestepManagerIterative& Picard = MDTM.dtManagerI(0);
    timestepManagerIterative& Newton = MDTM.dtManagerI(1);
    hEqn.checkSteadyConfig(Picard.tolerance(), Newton.tolerance());

    while (runTime.run())
    {
        if (!steady)
        {
            if (sourceEvent->isPresent()) sourceEvent->updateIndex(runTime.timeOutputValue());
            forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
            MDTM.updateDt();
        }

        runTime++;

noConvergence :
        Info << "Time = " << runTime.timeName() << nl << endl;

        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);

        if (sourceEvent->isPresent())
        {
            sourceEvent->updateValue(runTime);
            pmModel->sourceTerm() = sourceEvent->dtValuesAsField();
        }
        hEqn.updateSeepage();

        hEqn.solvePmModel(flowModels::RichardsEqn::couplingStep::BEFORE);

        //--- 1) Picard loop
        bool converged = Picard.solveEquation(hEqn, steady, 0);
        if (!converged)
        {
            hEqn.noConvergence(MDTM, runTime, 0);
            goto noConvergence;
        }

        //--- 2) Newton loop
        converged = Newton.solveEquation(hEqn, steady, 1);
        if (!converged)
        {
            hEqn.noConvergence(MDTM, runTime, 1);
            goto noConvergence;
        }

        hEqn.solvePmModel(flowModels::RichardsEqn::couplingStep::AFTER);

        //--- Compute variations
        if (!steady) MDTM.updateAllDerivatives();

        hEqn.info();

        if (steady)
        {
            runTime.write();
            if (writeResiduals)
                if (Pstream::master()) residualFile << runTime.timeName() << " " << MDTM.residual() << endl;
            if (Picard.residualDelta().first() < Picard.tolerance() && Newton.residualDelta().first() < Newton.tolerance())
            {
                runTime.writeAndEnd();
            }
        }
        else
        {
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
