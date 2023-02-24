.. _groundwater2DFoam:

groundwater2DFoam solver
========================

Description
-----------

Transient (or steady) solver for free-surface flow though porous media (2D modeling)

Time step is managed using an estimation of the time truncation error (see below).

Solver requires several geographical information (see below).

Infiltration over the domain can be constant and uniform (*infiltration* specified in **constant/transportProperties**), constant and non uniform (*infiltration* as a volScalarField in **0/infiltration** or variable and uniform (using an :ref:`infiltrationEventFile`)

A seepage option can be activated with Digital Elevation Model to constrain maximal water height (if potential > potentialDEM) and cell flux < 0 : solver sets potential=potentialDEM at this cell and computes exfiltration outflow)

The file **constant/potentialDEM** is a *volScalarField* which gives the maximal potential values used for seepage. This file can be generated automatically by the solver using a **.xyz** ASCII DEM file (entry *fileDEM* in **constant/transportProperties**).

Explicit fixed potentials (locations and values) can be given as a list in **constant/transportProperties**. The user-defined potential values in list can be overriden using the field **constant/potentialDEM** (requires activation of *useDEMtoFixPotential* option in **constant/transportProperties**).

To avoid numerical issues, a parameter *hWaterMin* impose the minimal water height for each cell. Its default is 0.01 and can be changed in **constant/transportProperties**. 

Configuration files
-------------------

**constant/transportProperties:**

.. code::

    // fluid and medium properties
    eps eps [0 0 0 0 0 0 0] 0.27; // porosity (can be volScalarField in constant/)
    rho rho [1 -3 0 0 0 0 0] 1000; // density
    mu mu [1 -1 -1 0 0 0 0] 1e-03; // viscosity

    hWaterMin [0 1 0 0 0 0 0] 0.01; // minimal water height for 'dry' cells

    infiltration infiltration [0 1 -1 0 0 0 0] -5e-9; // constant uniform infiltration velocity (optional)

    infiltrationEventFile "infiltrationOverTime.dat"; // to specify time-varying infiltration (uniform or non-uniform)

    seepage yes; // Seepage activation (requires constant/potentialDEM field)

    fileDEM dem.xyz; //(optional) mnt file read by solver to construct constant/potentialDEM if not present

    fixedPotentialList // list of the fixed potential positions (and user-defined potential values)
    (
        ( (584256.50 2413589.00 1) 145.2)
        ( (584556.50 2413889.00 1) 149.3)
    );

    useDEMtoFixPotential yes; // (use potentialDEM values to impose potential (instead of user-defined potential values)

**system/fvSolution:**

.. code::

  solvers
  {
      potential
      {
          solver          PCG;
          preconditioner  DIC;
          tolerance       1e-12;
          relTol          0;
      }

  }

**system/controlDict:**

.. code::

    adjustTimeStep yes;

    truncationError 0.001; // Allowed time-scheme truncation error used to manage time-step

    CSVoutput       true; // active the CmassBalance.csv output

    eventTimeTracking false; // to force the solver to compute solutions at each event time (patch/source/output)

    outputEventFile "output_times.dat"; // to specify output times (replace usual OpenFOAM writer based on writeInterval)

Required fields
---------------

- **0/potential :** The potential field

- **0/U :** The velocity field

- **constant/K :** permeability field

- **constant/z0 :** height of the impermeable layer

Optional fields
---------------

- **0/eps :** porosity field (replace eps in transport properties)

- **0/infiltration :** infiltration velocity field (replace infiltration in transport properties)

- **constant/potentialDEM :** field used by forcing potential modes: *fixedPotential* and *seepage*

Output computed fields
----------------------

- **0/hwater :** water height in the watershed

Timestep managing
-----------------

The computation of timestep for next iteration is directly computed using a user-defined *truncationError* related to the time scheme defined (**Euler**, **backward**, **CrankNicolson**). The time step formula for **backward** time-scheme is for example :

.. code::

  deltaT = Foam::pow(3 x truncationError x Hmax/dH3dT3max,1./3.)

where **dH3dT3max** is the maximal value of the thid time derivative and **Hmax** the value of hwater in this cell.

Steady simulation
-----------------

Solver can be run in *steady* mode using **-steady** option. The residual control and the under-relaxation factor should be specificied in **system/fvSolution** as :

.. code::

    residualControl
    {
        potential          1e-8;
    }

    relaxationFactors
    {
        equations
        {
            potential          0.5;
        }
    }

Seepage can be activated in *steady* mode.

