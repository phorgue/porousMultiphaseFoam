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
    groundwaterTransportFoam

Description
    Transient solver for Richards equation coupled with scalar transport 
    A Picard loop is used for linearization.
    Permeability is isotropic (K == volScalarField)

\*---------------------------------------------------------------------------*/

#include "PMFversion.H"
#include "fvCFD.H"
#include "multiMesh.H"
#include "dynamicRefineFvMesh.H"
#include "incompressiblePhase.H"
#include "twophasePorousMediumModel.H"
#include "porousMediumTransportModel.H"
#include "capillarityModel.H"
#include "relativePermeabilityModel.H"
#include "multiscalarMixture.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "patchEventFile.H"
#include "eventInfiltration.H"
#include "eventFlux.H"
#include "multiDtManager.H"
#include "RichardsEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    Foam::argList args(argc, argv);
    bool steady = false;
    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }
    PMFversion solverV;

    Info << "Create time\n" << Foam::endl;
    Time runTime(Time::controlDictName, args);

    //- Create mesh (simple or dual)
    autoPtr<dynamicFvMesh> meshPtrFluid(dynamicFvMesh::New(args, runTime));
    dynamicFvMesh& mesh = meshPtrFluid.ref();

    Info<< "Reading transportProperties" << endl;
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
    eventFlux::setEventFileRegistry(&patchEventList, "C");

    autoPtr<multiMesh> mMeshPtr(multiMesh::New(runTime, mesh, transportProperties));
    dynamicFvMesh& meshT = mMeshPtr.ref().fineMesh();

    //- fluid phase model
    autoPtr<incompressiblePhase> fluidPhase = incompressiblePhase::New(mesh, transportProperties, "theta");
    //- two-phase flow model
    autoPtr<twophasePorousMediumModel> pmModel =
            twophasePorousMediumModel::New("theta", mesh, transportProperties, fluidPhase);
    //- Richards' equation
    flowModels::RichardsEqn hEqn(mesh, transportProperties, pmModel(), fluidPhase(), steady);
    //- Transport's equation
    Info<< "Reading composition" << endl;
    autoPtr<porousMediumTransportModel> pmTransportModel = porousMediumTransportModel::New("theta", meshT, transportProperties);
    multiscalarMixture& composition = pmTransportModel->composition();

    //- Add
    volScalarField& thetaT = mMeshPtr->addField(hEqn.theta());
    volVectorField& UthetaT = mMeshPtr->addField(fluidPhase->U());
    surfaceScalarField& phiT = mMeshPtr->addField(hEqn.phi());

    //- create source events
    autoPtr<sourceEventFile> waterSourceEvent = sourceEventFile::New("sourceEventFileWater", transportProperties);
    waterSourceEvent->init(runTime,hEqn.h().name(), mesh, pmModel->sourceTerm().dimensions());
    List<sourceEventFile*>& tracerSourceEventList = pmTransportModel->sourceEventList();
    forAll(tracerSourceEventList,sourceEventi) tracerSourceEventList[sourceEventi]->init(runTime);
    forAll(patchEventList,patchEventi) patchEventList[patchEventi]->init(runTime);

    //- create time managers
    multiDtManager MDTM(runTime, tracerSourceEventList, patchEventList);
    MDTM.addIterativeAlgorithm(hEqn.theta(), "Picard");
    MDTM.addIterativeAlgorithm(hEqn.theta(), "Newton");
    const labelList* fixedPotentialIDListPtr  = &hEqn.seepageIDList();
    MDTM.addField(hEqn.h(), fixedPotentialIDListPtr);
    forAll(composition.Y(), speciesi) MDTM.addField(composition.Y()[speciesi]);
    if (meshT.dynamic()) MDTM.setDynamicMesh(true);

    //-Output event
    autoPtr<outputEventFile> outputEventF = outputEventFile::New(runTime, mesh);
    outputEventF->addField(hEqn.h(), hEqn.phi());
    outputEventF->addField(hEqn.theta(), hEqn.phi(), "waterMassBalance.csv", true);
    autoPtr<outputEventFile> outputEventT = outputEventFile::New(runTime, meshT);
    forAll(composition.Y(), speciei) {
        outputEventT->addField(composition.Y()[speciei], phiT, thetaT, composition.R(speciei), composition.Y()[speciei].name()+"massBalance.csv");
    }
    outputEventF->init();
    outputEventT->init();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    timestepManagerIterative& Picard = MDTM.dtManagerI(0);
    timestepManagerIterative& Newton = MDTM.dtManagerI(1);

    while (runTime.run())
    {
        if (waterSourceEvent->isPresent())  waterSourceEvent->updateIndex(runTime.timeOutputValue());
        forAll(tracerSourceEventList,tracerSourceEventi) tracerSourceEventList[tracerSourceEventi]->updateIndex(runTime.timeOutputValue());
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());

        MDTM.updateDt();

        runTime++;

noConvergence :
        Info << "Time = " << runTime.timeName() << nl << endl;

        //- Update source term
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
        if (waterSourceEvent->isPresent())
        {
            waterSourceEvent->updateValue(runTime);
            pmModel->sourceTerm() = waterSourceEvent->dtValuesAsField();
        }
        hEqn.updateSeepage();

        Tuple2<scalar, scalar> residualDelta(1.00001, 0);

        //- 1) Richard's equation (Picard loop)
        Picard.reset();
        while (residualDelta.first() > Picard.tolerance() && Picard.iter() != Picard.maxIter())
        {
            Picard++;
            Info << "*** Picard iteration " << Picard.iter() << endl;
            residualDelta = hEqn.solvePicard(Picard.tolerance());
            hEqn.updateProperties(false);
            if (residualDelta.first() > 10)
            {
                Warning() << "Non-physical values reached, reducing time step by factor dTFactDecrease" << nl << endl;
                hEqn.noConvergence(MDTM, runTime, 0);
                goto noConvergence;
            }
        }
        if (residualDelta.first()  > Picard.tolerance())
        {
            Info << endl;
            if (MDTM.adjustTimeStep()) Warning() << " Max iteration reached in Picard loop, reducing time step by factor dTFactDecrease" << nl << endl;
            else FatalErrorIn("groundwaterTransportFoam.C") << "Non-convergence of Picard algorithm with fixed timestep => Decrease the time step or increase tolerance" << exit(FatalError);
            hEqn.noConvergence(MDTM, runTime, 0);
            goto noConvergence;
        }

        //--- 2) Newton loop
        Newton.reset();
        while (residualDelta.first() > Newton.tolerance() && Newton.iter() != Newton.maxIter())
        {
            Newton++;
            Info << "*** Newton iteration " << Newton.iter() << endl;
            residualDelta = hEqn.solveNewton(Newton.tolerance());
            hEqn.updateProperties(true);
            if (residualDelta.first() > 10)
            {
                Warning() << "Non-physical values reached, reducing time step by factor dTFactDecrease" << nl << endl;
                hEqn.noConvergence(MDTM, runTime, 1);
                goto noConvergence;
            }
        }
        if (residualDelta.first() > Newton.tolerance())
        {
            Info << endl;
            if (MDTM.adjustTimeStep()) Warning() <<  " Max iteration reached in Newton loop, reducing time step by factor dTFactDecrease" << nl << endl;
            else FatalErrorIn("groundwaterFoam.C") << "Non-convergence of Newton algorithm with fixed timestep => Decrease the time step or increase tolerance" << exit(FatalError);
            hEqn.noConvergence(MDTM, runTime, 1);
            goto noConvergence;
        }
        hEqn.info();

        //- 3) scalar transport
        forAll(tracerSourceEventList,tracerSourceEventi) tracerSourceEventList[tracerSourceEventi]->updateValue(runTime);
        if (mMeshPtr->dynamic()) composition.updateNormalizedGradY();
        bool hasChanged = mMeshPtr->update();
        if (hasChanged) forAll(tracerSourceEventList,tracerSourceEventi) tracerSourceEventList[tracerSourceEventi]->onMeshChanged();
        pmTransportModel->solveTransport(UthetaT, phiT, thetaT, pmModel->exchangeTerm());

        //- C and water mass balance computation
        MDTM.updateAllDerivatives();
        outputEventF->write();
        outputEventT->write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
