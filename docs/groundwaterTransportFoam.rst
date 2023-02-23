.. _groundwaterTransportFoam:

groundwaterTransportFoam solver
===============================

Description
-----------

This solver is a merge between :ref:`groundwaterFoam` and :ref:`porousScalarTransportFoam` (see dedicated wiki pages for details)

At each time iteration, it solves

    - the Richards' equation for the water flow
    - passive scalar transport equation for the species

Configuration files
-------------------

**constant/transportProperties :**

.. code::

    Ss Ss [0 -1 0 0 0 0 0] 0; // specific storage

    phase.theta
    {
        rho rho [1 -3 0 0 0 0 0] 1e3; // density of the fluid
        mu mu [1 -1 -1 0 0 0 0] 1e-3; // dynamic viscosity of the fluid
    }

    relativePermeabilityModel  VanGenuchten;

    capillarityModel	VanGenuchten;

    VanGenuchtenCoeffs // parameters of the Van Genuchten kr/pc model
    {
        thetamin 0.102; // minimal saturation
        thetamax 0.368; // maximal saturation
        m 0.5;
        alpha 3.35;
    }

    eps eps [0 0 0 0 0 0 0] 0.25; // porosity (can be volScalarField in constant/)

    Dm Dm [0 2 -1 0 0 0 0] 1e-9; // molecular diffusion

    porousTransport // specific dictionary for transport
    {
        phaseName a; // to specify the flux field (phia here) and velocity field (Ua)
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

    eventFileTracerSource injection.dat; // to specify event file for tracer source term
    eventFileWaterSource water_injection.dat; // to specify event file for water source term

**system/controlDict :**

.. code::

    adjustTimeStep yes;

    timeStepControl Picard; // Picard or dthetamax

    //- for h variation time step control
    dthetamax           0.005; // theta variation instruction for computing the time step

    //- for time step control
    truncationError 0.001; // global truncation error used to manage time-step
    truncationError_C 0.001; // (optional) tracer only truncation error
    truncationError_h 0.001; // (optional) potential only truncation error

    CSVoutput       true; // active the waterMassBalance.csv and CmassBalance.csv outputs

    eventTimeTracking false; // to force the solver to compute solutions at each event time (patch/source/output)


Required fields
---------------

- **0/h :** The potential field
- **0/Utheta :** The velocity field
- **0/C :** The concentration field
- **constant/g :** gravity field
- **constant/K :** permeability field

Optional fields
---------------

- Other spatially defined parameters : **alphaL** , **alphaT** , **eps**.
- **0/thetamin** and **0/thetamax :** spatialized minimal and maximal saturation (replaces thetamin and thetamax in **transportProperties**)

Timestep managing
-----------------

The timestep is managed as for the two original solvers, taking the minimal deltaT required by water transport and scalar transport (with eventually different truncation error parameters for **C** and **h**)

See :ref:`groundwaterFoam` and :ref:`porousScalarTransportFoam` for more information.
