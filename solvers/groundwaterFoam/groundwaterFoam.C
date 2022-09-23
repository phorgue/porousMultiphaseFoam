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
#include "porousMediumModel.H"
#include "capillarityModel.H"
#include "relativePermeabilityModel.H"
#include "fixedValueFvPatchField.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "patchEventFile.H"
#include "eventInfiltration.H"
#include "timestepManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addBoolOption("steady", "to run steady flow simulation");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "readPicardNewtonControls.H"
    #include "readTimeControls.H"
    #include "createthetaFields.H"
    #include "readEvent.H"
    #include "readForcing.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    bool steady = args.found("steady");
    if (steady) maxIterPicard = 1;
    label iterPicard=0;
    label iterNewton=0;

    while (runTime.run())
    {
        if (!steady)
        {
            if (sourceEventIsPresent) sourceEvent.updateIndex(runTime.timeOutputValue());
            forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
            #include "setDeltaT.H"
        }

        runTime++;

noConvergence :
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "computeSourceTerm.H"
        scalar deltahIter = 1;
        scalar hEqnResidual = 1.00001;
        scalar hEqnResidualSigned = 0;

        //--- 1) Picard loop
        iterPicard = 0;
        while ( hEqnResidual > tolerancePicard && iterPicard != maxIterPicard )
        {
            iterPicard++;
            #include "hEqnPicard.H"
            #include "checkResidual.H"
            Info << "Picard iteration " << iterPicard << ": max(deltah) = " << deltahIter << ", residual = " << hEqnResidualSigned << endl;
        }
        if ( !steady &&  hEqnResidual > tolerancePicard )
        {
            Info << endl;
            if (adjustTimeStep) Warning() << " Max iteration reached in Picard loop, reducing time step by factor dTFactDecrease" << nl << endl;
            else FatalErrorIn("groundwaterFoam.C") << "Non-convergence of Picard algorithm with fixed timestep => Decrease the time step or increase tolerance" << exit(FatalError);
            #include "rewindTime.H"
            goto noConvergence;
        }

        //--- 2) Newton loop
        iterNewton = 0;
        while ( hEqnResidual > toleranceNewton && iterNewton != maxIterNewton)
        {
            iterNewton++;
            #include "hEqnNewton.H"
            #include "checkResidual.H"
            Info << "Newton iteration : " << iterNewton << ": max(deltah) = " << deltahIter << ", residual = " << hEqnResidualSigned << endl;
        }
        if ( !steady && hEqnResidual > toleranceNewton )
        {
            Info << endl;
            if (adjustTimeStep) Warning() <<  " Max iteration reached in Newton loop, reducing time step by factor dTFactDecrease" << nl << endl;
            else FatalErrorIn("groundwaterFoam.C") << "Non-convergence of Newton algorithm with fixed timestep => Decrease the time step or increase tolerance" << exit(FatalError);
            #include "rewindTime.H"
            goto noConvergence;
        }

        //--- Compute variations
        if (!steady) dtManager.updateDerivatives();
        scalarField dtheta_tmp = mag(theta.internalField()-theta.oldTime().internalField());
        scalar dtheta = gMax(dtheta_tmp);

        Info << "Saturation theta: " << " Min(theta) = " << gMin(theta.internalField()) << " Max(theta) = " << gMax(theta.internalField()) << " dthetamax = " << dtheta << endl;
        Info << "Head pressure h: " << " Min(h) = " << gMin(h.internalField()) << " Max(h) = " << gMax(h.internalField()) << endl;

        #include "waterMassBalance.H"
        if (steady)
        {
            if (hEqnResidual < tolerancePicard) runTime.writeAndEnd();
        }
        else
        {
            #include "eventWrite.H"
        }
 
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
