.. _groundwaterTransport2DFoam:

groundwaterTransport2DFoam solver
=================================


Description
-----------

This solver is a merge between :ref:`groundwater2DFoam` and :ref:`porousScalarTransport2DFoam` (see dedicated wiki pages for details).

At each time iteration, it solves

    - the free-surface flow equation in porous media
    - passive scalar transport equation for the species

Configuration files
-------------------

**constant/transportProperties :**

.. code::

    // fluid and medium properties
    eps eps [0 0 0 0 0 0 0] 0.27; // porosity (can be volScalarField in constant/)
    rho rho [1 -3 0 0 0 0 0] 1000; // density
    mu mu [1 -1 -1 0 0 0 0] 1e-03; // viscosity

    infiltration infiltration [0 1 -1 0 0 0 0] -5e-9; // constant uniform infiltration velocity

    eventFileInfiltration "infiltrationOverTime.dat"; // to specify time-varying infiltration (uniform or non-uniform)

    seepage yes; // Seepage activation (requires constant/potentialDEM field)

    fileDEM dem.xyz; //(optional) mnt file read by solver to construct constant/potentialDEM if not present

    fixedPotentialList // list of the fixed potential positions (and user-defined potential values)
    (
        ( (584256.50 2413589.00 1) 145.2)
        ( (584556.50 2413889.00 1) 149.3)
    );

    useDEMtoFixPotential yes; // (use potentialDEM values to impose potential (instead of user-defined potential values)

    Dm Dm [0 2 -1 0 0 0 0] 1e-9; // molecular diffusion

    porousTransport // specific dictionary for transport
    {
        Kd Kd [-1 3 0 0 0 0 0] 1e-3;
        rs rs [1 -3 0 0 0 0 0] 1000;
        epsTotal epsTotal [0 0 0 0 0 0 0] 0.30;
        lambda lambda [0 0 -1 0 0 0 0 ] 0;// decay of the C scalar
    }

    dispersionModel alphaDispersion; // dispersion model

    alphaDispersionCoeffs
    {
        tau tau [0 0 0 0 0 0 0] 2; // tortuosity
        alphaL alphaL [0 1 0 0 0 0 0] 0.01; // longitudinal dispersivity
        alphaT alphaT [0 1 0 0 0 0 0] 0.002; // transverse dispersivity
    }

    eventFileTracerSource injection.dat; // to specify event file for time-dependent source term

**system/controlDict :**

.. code::

    adjustTimeStep yes;

    //- for Scalar transport time step control
    truncationError 0.001; // global truncation error used to manage time-step
    truncationError_potential 0.001; // (optional) potential only truncation error used to manage time-step
    truncationError_C 0.001; // (optional) tracer only truncation error used to manage time-step

    CSVoutput       true; // active the CmassBalance.csv output

    eventTimeTracking true; // to force the solver to compute solutions at each event time (patch/source/output)


Required fields
---------------

- **0/potential :** The potential field
- **0/U :** The velocity field
- **0/C :** The concentration field
- **constant/z0 :** height of the impermeable layer
- **constant/K :** permeability field

Optional fields
---------------

- Other spatially defined parameters : **alphaL** , **alphaT** , **eps**.
- **0/infiltration :** infiltration velocity field (replace infiltration in transport properties)

Timestep managing
-----------------

The timestep is managed as for the two original solvers, taking the minimal deltaT required by water transport and scalar transport (with eventually different truncation error parameters for **C** and **potential**).

See :ref:`groundwater2DFoam` and :ref:`porousScalarTransport2DFoam` for more information.
