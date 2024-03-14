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

#include "fvCFD.H"
#include "harmonic.H"
#include "incompressiblePhase.H"
#include "twophasePorousMediumModel.H"
#include "capillarityModel.H"
#include "relativePermeabilityModel.H"
#include "fixedValueFvPatchField.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "patchEventFile.H"
#include "eventInfiltration.H"
#include "multiDtManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addBoolOption("steady", "to run steady flow simulation");

    Foam::argList args(argc, argv);
    bool steady = args.found("steady");

    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }
    #include "../headerPMF.H"

    Info << "Create time\n" << Foam::endl;
    Time runTime(Time::controlDictName, args);

    #include "createMesh.H"

    Info<< "\nReading g" << endl;
    const meshObjects::gravity& g = meshObjects::gravity::New(runTime);

    #include "createFields.H"
    bool massConservative = transportProperties.lookupOrDefault<bool>("massConservative",true);
    #include "readForcing.H"

    #include "createthetaFields.H"

    //- create source event for water 
    autoPtr<sourceEventFile> sourceEvent = sourceEventFile::New("sourceEventFileWater", transportProperties);
    sourceEvent->init(runTime, h.name(), mesh, sourceTerm.dimensions());

    //- create time managers
    labelList* fixedPotentialIDListPtr  = &fixedPotentialIDList;
    List<sourceEventFile*> sourceEventList(1, sourceEvent);
    multiDtManager MDTM(runTime, sourceEventList, patchEventList);
    MDTM.addIterativeAlgorithm(theta, "Picard", steady);
    MDTM.addIterativeAlgorithm(theta, "Newton", steady);
    MDTM.addField(h, fixedPotentialIDListPtr);

    //- output event
    autoPtr<outputEventFile> outputEvent = outputEventFile::New(runTime, mesh);
    outputEvent->addField(h, phi);
    outputEvent->addField(theta, phi, "waterMassBalance.csv", true);



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    timestepManagerIterative& Picard = MDTM.dtManagerI(0);
    timestepManagerIterative& Newton = MDTM.dtManagerI(1);

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
            sourceTerm = sourceEvent->dtValuesAsField();
        }
        #include "updateForcing.H"

        scalar deltahIter = 1;
        scalar hEqnResidualN = 1.00001;
        scalar hEqnResidualP = 1.00001;
        scalar hEqnResidualInit = 1.00001;

        //--- 1) Picard loop
        Picard.reset();
        while ( hEqnResidualInit > Picard.tolerance() && Picard.iter() != Picard.maxIter() )
        {
            Picard++;
            #include "hEqnPicard.H"
            #include "updateProperties.H"
            #include "computeResidualN.H"
            Info << "Picard iteration " << Picard.iter() << ": max(deltah) = " << deltahIter << ", residualP = " << hEqnResidualP << ", residualN = " << hEqnResidualN << endl;
            if ( hEqnResidualInit > 10)
            {
                Warning() << "Non-physical values reached, reducing time step by factor dTFactDecrease" << nl << endl;
                Picard.reset(Picard.maxIter());
                #include "rewindTime.H"
                goto noConvergence;
            }
        }
        if ( !steady &&  hEqnResidualInit > Picard.tolerance() )
        {
            if (MDTM.adjustTimeStep()) Warning() << " Max iteration reached in Picard loop, reducing time step by factor dTFactDecrease" << nl << endl;
            else FatalErrorIn("groundwaterFoam.C") << "Non-convergence of Picard algorithm with fixed timestep => Decrease the time step / add field h relaxation / increase number of Picard iterations" << exit(FatalError);
            #include "rewindTime.H"
            goto noConvergence;
        }

        //--- 2) Newton loop
        Newton.reset();
        while ( hEqnResidualN > Newton.tolerance() && Newton.iter() != Newton.maxIter())
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
            Info << "Newton iteration : " << Newton.iter() << ": max(deltah) = " << deltahIter << ", residualN = " << hEqnResidualN << endl;
            if ( hEqnResidualN > 10)
            {
                Warning() << "Non-physical values reached, reducing time step by factor dTFactDecrease" << nl << endl;
                Newton.reset(Newton.maxIter());
                #include "rewindTime.H"
                goto noConvergence;
            }
        }
        if ( !steady && hEqnResidualN > Newton.tolerance() )
        {
            Info << endl;
            if (MDTM.adjustTimeStep()) Warning() <<  " Max iteration reached in Newton loop, reducing time step by factor dTFactDecrease" << nl << endl;
            else FatalErrorIn("groundwaterFoam.C") << "Non-convergence of Newton algorithm with fixed timestep => Decrease the time step or increase tolerance" << exit(FatalError);
            #include "rewindTime.H"
            goto noConvergence;
        }

        //--- Compute variations
        if (!steady) MDTM.updateAllDerivatives();
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

        if (steady)
        {
            runTime.write();
            if (writeResiduals)
            {
                OFstream residualFile("residuals.csv", IOstreamOption(), true);
                residualFile << runTime.timeName() << " " << mag(hEqnResidualP) << endl;
            }
            if (hEqnResidualInit < Picard.tolerance()) runTime.writeAndEnd();
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
