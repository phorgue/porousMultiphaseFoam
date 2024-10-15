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
    Transient solver for Richards equation coupled with scalar transport 
    A Picard loop is used for linearization.
    Permeability is isotropic (K == volScalarField)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "multiMesh.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
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
#include "processorFvPatchField.H"
#include "symmetryPlanePolyPatch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addBoolOption("dualDynamicMesh", "to run steady flow simulation");

    Foam::argList args(argc, argv);
    bool steady = false;
    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }

    #include "../headerPMF.H"
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
    autoPtr<multiMesh> mMeshPtr(multiMesh::New(mesh, transportProperties));
    dynamicFvMesh& meshT = mMeshPtr.ref().fineMesh();

    Info << "\nReading g" << endl;
    const meshObjects::gravity& g = meshObjects::gravity::New(runTime);

    #include "createFields.H"
    volScalarField& thetaT = mMeshPtr->addField(theta);
    volVectorField& UthetaT = mMeshPtr->addField(Utheta);
    surfaceScalarField& phiT = mMeshPtr->addField(phi);

    bool massConservative = transportProperties.lookupOrDefault<bool>("massConservative",true);
    #include "readForcing.H"

    bool writeResiduals = false; //- stationary run not possible with transport
    #include "createthetaFields.H"

    autoPtr<sourceEventFile> waterSourceEvent = sourceEventFile::New("sourceEventFileWater", transportProperties);
    waterSourceEvent->init(runTime, h.name(), mesh, sourceTerm.dimensions());
    forAll(tracerSourceEventList,sourceEventi) tracerSourceEventList[sourceEventi]->init(runTime);
    forAll(patchEventList,patchEventi) patchEventList[patchEventi]->init(runTime);

    //- create time managers
    multiDtManager MDTM(runTime, tracerSourceEventList, patchEventList);
    MDTM.addIterativeAlgorithm(theta, "Picard");
    MDTM.addIterativeAlgorithm(theta, "Newton");
    forAll(composition.Y(), speciesi) MDTM.addField(composition.Y()[speciesi]);
    if (meshT.dynamic()) MDTM.setDynamicMesh(true);

    //-Output event
    autoPtr<outputEventFile> outputEventF = outputEventFile::New(runTime, mesh);
    outputEventF->addField(h, phi);
    outputEventF->addField(theta, phi, "waterMassBalance.csv", true);
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
        if (waterSourceEvent->isPresent())
        {
            waterSourceEvent->updateValue(runTime);
            sourceTerm = waterSourceEvent->dtValuesAsField();
        }
        #include "updateForcing.H"

        scalar deltahIter = 1;
        scalar hEqnResidualMax = 1.00001;
        scalar hEqnResidualInit = 1.00001;

        //- 1) Richard's equation (Picard loop)
        Picard.reset();
        while ( hEqnResidualInit > Picard.tolerance() && Picard.iter() != Picard.maxIter() )
        {
            Picard++;
            #include "hEqnPicard.H"
            #include "updateProperties.H"
            #include "computeResidualN.H"
            Info << "Picard iteration " << Picard.iter() << ": max(deltah) = " << deltahIter << ", max(residual) = " << hEqnResidualMax << endl;
            if ( hEqnResidualInit > 10)
            {
                Warning() << "Non-physical values reached, reducing time step by factor dTFactDecrease" << nl << endl;
                Picard.reset(Picard.maxIter());
                #include "rewindTime.H"
                goto noConvergence;
            }
        }
        if ( hEqnResidualInit > Picard.tolerance() )
        {
            Info << endl;
            if (MDTM.adjustTimeStep()) Warning() << " Max iteration reached in Picard loop, reducing time step by factor dTFactDecrease" << nl << endl;
            else FatalErrorIn("groundwaterTransportFoam.C") << "Non-convergence of Picard algorithm with fixed timestep => Decrease the time step or increase tolerance" << exit(FatalError);
            #include "rewindTime.H"
            goto noConvergence;
        }

        //--- 2) Newton loop
        Newton.reset();
        while ( hEqnResidualMax > Newton.tolerance() && Newton.iter() != Newton.maxIter())
        {
            if (Picard.iter() == 0)
            {
                #include "computeResidualN.H"
                Picard++;
            }
            Newton++;
            #include "hEqnNewton.H"
            #include "updateProperties.H"
            #include "computeResidualN.H"
            Info << "Newton iteration : " << Newton.iter() << ": max(deltah) = " << deltahIter << ", max(residual) = " << hEqnResidualMax << endl;
            if ( hEqnResidualMax > 10)
            {
                Warning() << "Non-physical values reached, reducing time step by factor dTFactDecrease" << nl << endl;
                Newton.reset(Newton.maxIter());
                #include "rewindTime.H"
                goto noConvergence;
            }
        }
        if ( !steady && hEqnResidualMax > Newton.tolerance() )
        {
            Info << endl;
            if (MDTM.adjustTimeStep()) Warning() <<  " Max iteration reached in Newton loop, reducing time step by factor dTFactDecrease" << nl << endl;
            else FatalErrorIn("groundwaterFoam.C") << "Non-convergence of Newton algorithm with fixed timestep => Decrease the time step or increase tolerance" << exit(FatalError);
            #include "rewindTime.H"
            goto noConvergence;
        }

        //--- Compute variations
        scalarField dtheta_tmp = mag(theta.internalField()-theta.oldTime().internalField());
        scalar dtheta = gMax(dtheta_tmp);

        //- water mass balance terminal display
        Info << "Saturation theta: " << " Min(theta) = " << gMin(theta.internalField()) << " Max(theta) = " << gMax(theta.internalField()) << " dthetamax = " << dtheta << endl;
        Info << "Head pressure h: " << " Min(h) = " << gMin(h.internalField()) << " Max(h) = " << gMax(h.internalField()) << endl;
        Info << "Water mass balance (m3/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value() << " ; ";
        forAll(phi.boundaryField(),patchi)
        {
            if (mesh.boundaryMesh()[patchi].type() == "patch")
            {
                Info << phi.boundaryField()[patchi].patch().name() << " = " <<  gSum(phi.boundaryField()[patchi]) << " ; ";
            }
        }
        Info << endl;

        //- 3) scalar transport
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
        forAll(tracerSourceEventList,tracerSourceEventi) tracerSourceEventList[tracerSourceEventi]->updateValue(runTime);

        if (mMeshPtr->dynamic()) composition.updateNormalizedGradY();
        mMeshPtr->update();
        if (meshT.changing()) forAll(tracerSourceEventList,tracerSourceEventi) tracerSourceEventList[tracerSourceEventi]->onMeshChanged();
        pmTransportModel->solveTransport(UthetaT, phiT, thetaT, porousModel->exchangeTerm());

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
