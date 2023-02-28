.. _groundwaterFoam:

groundwaterFoam solver
======================

Description
-----------

Transient (or steady) solver for groundwater flow in porous media (Richards' equation with isotropic permeability **volScalarField**)

The non-linearities of relative permeability and capillary pressure are handled at each timestep by iterative methods:

1) Picard's iterations first to converge under a user-defined *Picard tolerance*
2) Newton's iterations at each time step to converge under a user-defined *Newton tolerance*

Each algorithm can be skipped by setting *tolerance* equal to 1. If the convergence is not reached in **maxIter** iterations, the current time iteration starts again with a smaller timestep.

Time step is managed using an estimation of the time truncation error of the saturation variation (see below).

Relative permeabilites and capillary model can be chose in the :ref:`porousModels`.

Seepage can be activated by specifying the name of the patch *patchDEM* (representing the top of the watershed). At given time iteration, if a neighbour cell is saturated and cell flux < 0 : solver sets h=dz/2 at this cell and computes exfiltration outflow.

.. warning::
   The Newton's algorithm does not handle correctly the seepage function and should be deactivated when seepage occurs in a simulation (by setting tolerance Newton equal to 1)

This solver allows the use of *dualPorosity* model, see :ref:`porousMediumModels` for details.

Configuration files
-------------------

**constant/transportProperties :**

.. code::

    Ss Ss [0 -1 0 0 0 0 0] 0; // specific storage

    patchDEM nameOfThePatch; // patch used for seepage options

    phase.theta
    {
        rho	rho [1 -3 0 0 0 0 0] 1e3; // density of the fluid
        mu	mu [1 -1 -1 0 0 0 0] 1e-3; // dynamic viscosity of the fluid
    }

    relativePermeabilityModel  VanGenuchten;

    capillarityModel	VanGenuchten;

    VanGenuchtenCoeffs // parameters of the Van Genuchten kr/pc model
    {
        thetamin 0.102; // minimal saturation
        thetamax 0.368; // maximal saturation
        m	0.5;
        alpha 3.35;
    }

**system/fvSolution :**

.. code::

    solvers
    {
       "h" // for Picard's solver
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-12;
            relTol          0;
        }
        "deltah" // for Newton's solver
        {
            solver          PBiCG;
            preconditioner  DILU;
            tolerance       1e-12;
            relTol          0;
        }

    }

    Picard
    {
        tolerance 0.001; // tolerance for the Picard's algorithm
        maxIter 10; // maximum iterations for the Picard's algorithm
    }

    Newton
    {
        tolerance 1e-06; // tolerance for the Newton's algorithm
        maxIter 10; // maximum iterations for the Newton's algorithm
    }

**system/controlDict**

.. code::

    adjustTimeStep  yes;

    truncationError   0.005; // theta variation instruction for computing the time step

    CSVoutput       true; // active the waterMassBalance.csv output

    outputEventFile outputList.dat; // to specify the writing time outputs using linear time interpolations (replaces usual write() function of OpenFOAM)

    eventTimeTracking false; // to force the solver to  explicitly compute output event time solutions (instead of time linear interpolations)

Required fields
---------------

- **0/h :** The potential field
- **0/Utheta :** The velocity field
- **constant/g :** gravity field
- **constant/K :** permeability field

Optional fields
---------------

- **0/thetamin** and **0/thetamax :** spatialized minimal and maximal saturation (replace *thetamin* and *thetamax* in **transportProperties**)

- **0/m** and **0/alpha :** spatialized Van Genuchten parameters (replace *m* and *alpha* in **transportProperties**)

Timestep managing
-----------------

The computation of timestep for next iteration is directly computed using truncation error related to the time scheme used. Due to the Newton's  method time dicretization, only the *Euler* scheme is available with:

.. code::

  deltaT = Foam::pow(2 x truncationError x Hmax[speciesi]/dH2dT2max[speciesi],1./3.)

where **dH2dT2max** is the maximal value of the second order time derivative and **Hmax** the value of hwater in this cell.

Steady simulation
-----------------

Solver can be run in *steady* mode using **-steady** option. The under-relaxation factor on the pressure head should be specificied in **system/fvSolution** as :

.. code::

    relaxationFactors
    {
        fields
        {
            h          0.01;
        }
    }

Simulation occurs until *Picard's tolerance* is reached (note the Newton's algorithm is not used in *steady* mode).

Seepage can be activated in *steady* mode.
