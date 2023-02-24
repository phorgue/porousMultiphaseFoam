.. _porousScalarTransport2DFoam:

porousScalarTransport2DFoam solver
==================================

Description
-----------

Solves the transport equation of a passive scalar for free-surface flow in porous media with dispersion and retard coefficient model.

Multiple time-dependent injection source points can be specificied using the keyword *sourceEventFileTracer* in **constant/transportProperties** (see :ref:`sourceEventFile`).

The transport is implicitly solved and based on pre-computed flux **phi** or velocity **U** field and water depth **hwater**.

Time step is managed using an estimation of the time truncation error of the scalar variation.

Configuration files
-------------------

**constant/transportProperties :**

Example for one specie transport (or multi-species with similar properties):

.. code::

  eps eps [0 0 0 0 0 0 0] 0.25; // porosity (can be volScalarField in constant/)

  Dm Dm [0 2 -1 0 0 0 0] 1e-9; // molecular diffusion

  porousTransport // specific dictionary for transport
  {
    phaseName theta; // to specify the flux field (phitheta here) or velocity field (Utheta)
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

  sourceEventFileTracer injectionFile.dat;  // to specify event file for time-dependent source term

Example for multi-species properties can be found at :ref:`porousScalarTransportFoam` page.

**system/controlDict :**

.. code::

    adjustTimeStep yes;

    truncationError 0.001; // Allowed time-scheme truncation error used to manage time-step

    CSVoutput       true; // active the CmassBalance.csv output

    eventTimeTracking true; // to force the solver to compute solutions at each event time (patch/source/output)

Required fields
---------------

- **0/C :** The concentration field

- **0/hwater :** The pre-computed water depth

- **0/U :** Used to compute dispersion coefficient. Can be used to compute **0/phi** if not present (*Warning: velocity field U in FV formulation is usually reconstructed by interpolation from flux field phi and is not conservative, prefer the use of real phi field from flow field simulation*


Optional fields
---------------

- **0/phi:** The pre-computed flux field where `a` can be changed  in `porousTransport` dictionary (can be initialized using **0/U**)

- Other spatially defined parameters : **alphaL** , **alphaT**, **eps**

Timestep managing
-----------------

The computation of timestep for next iteration is directly computed using truncation error related to the time scheme defined (**Euler**, **backward**, **CrankNicolson**). The time step formula for **backward** time-scheme is for example :

.. code::

  deltaT = Foam::pow(3 x truncationError x Cmax[speciesi]/dC3dT3max[speciesi],1./3.)

where **dC3dT3maxmaximal** is the maximal value of the thid time derivative and **Cmax** the value of C in this cell._
