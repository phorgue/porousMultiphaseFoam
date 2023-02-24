.. _porousScalarTransportFoam:

porousScalarTransportFoam solver
================================

Description
-----------

Solves the transport equation for a passive scalar in porous media with dispersion and retard coefficient model.

Multiple time-dependent injection source points can be specificied using the keyword *sourceEventFileTracer* in **constant/transportProperties** (see :ref:`sourceEventFile`).

The transport is implicitly solved using a pre-computed flux **phi** or velocity **U** and optional water content **theta** or saturation **S**

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
    phaseName theta; // to specify the flux field (phitheta here) or velocity field (Utheta). U and phi is not precised.
    Kd Kd [-1 3 0 0 0 0 0] 1e-3;
    rs rs [1 -3 0 0 0 0 0] 1000;
    epsTotal epsTotal [0 0 0 0 0 0 0] 0.30;
    lambda lambda [0 0 -1 0 0 0 0 ] 1e-9;// decay of the C scalar
  }

  dispersionModel alphaDispersion; // dispersion model

  alphaDispersionCoeffs
  {
    tau tau [0 0 0 0 0 0 0] 2; // tortuosity
    alphaL alphaL [0 1 0 0 0 0 0] 0.01; // longitudinal dispersivity
    alphaT alphaT [0 1 0 0 0 0 0] 0.002; // transverse dispersivity
  }

  sourceEventFileTracer injection.dat; // to specify event file for time-dependent source term

Example for multi-species transport (C1/C2 species) :

.. code::

    eps eps [0 0 0 0 0 0 0]	0.25;

    porousTransport
    {
        phaseName oil;
    }

    species
    (
        C1
        C2
    );


    C1
    {
      Dm Dm [0 2 -1 0 0 0 0] 1e-9;
      porousTransport
      {
          Kd Kd [-1 3 0 0 0 0 0] 1e-3;
          rs rs [1 -3 0 0 0 0 0] 0;
          epsTotal epsTotal [0 0 0 0 0 0 0] 0.30;
          lambda lambda [0 0 -1 0 0 0 0 ] 1.1574e-6;
      }

      dispersionModel alphaDispersion;
      alphaDispersionCoeffs
      {
          tau tau [0 0 0 0 0 0 0] 2;
          alphaL alphaL [0 1 0 0 0 0 0] 0.01;
          alphaT alphaT [0 1 0 0 0 0 0] 0.002;
      }
    }

    C2
    {
      Dm Dm [0 2 -1 0 0 0 0] 1e-9;
      porousTransport
      {
          Kd Kd [-1 3 0 0 0 0 0] 1e-3;
          rs rs [1 -3 0 0 0 0 0] 0;
          epsTotal epsTotal [0 0 0 0 0 0 0] 0.30;
          lambda lambda [0 0 -1 0 0 0 0 ] 1.1574e-6;
      }

      dispersionModel alphaDispersion;
      alphaDispersionCoeffs
      {
          tau tau [0 0 0 0 0 0 0] 2;
          alphaL alphaL [0 1 0 0 0 0 0] 0.01;
          alphaT alphaT [0 1 0 0 0 0 0] 0.002;
      }
    }

**system/controlDict :**

.. code::

    adjustTimeStep yes;

    truncationError 0.001; // Allowed time-scheme truncation error used to manage time-step

    CSVoutput       true; // active the CmassBalance.csv output

    eventTimeTracking true; // to force the solver to compute solutions at each event time (patch/source/output)


Required fields
---------------

- **0/C :** The concentration field (or **0/C1** **0/C2** in the multi-specie example above)

- **0/UphaseName :**  (phaseName is read from **transportProperties**) Used to compute dispersion coefficient. Can be used to compute **0/phiphaseName** if not present (*Warning: velocity field U in FV formulation has been reconstructed from flux field phi and is not conservative, prefer the use of real phi field*


Optional fields
---------------

- **0/phiphaseName :** The pre-computed flux field where *phaseName* can be changed  in *porousTransport* dictionary (**0/phioil** in multi-specie example)

- **0/SphaseName :** Saturation field (S=1 if not present)

- **0/theta :** Water content (overwrites Saturation if both are present)

- Other spatially defined parameters : **alphaL** , **alphaT** , **eps**.

Timestep managing
-----------------

The computation of timestep for next iteration is directly computed using truncation error related to the time scheme defined (**Euler**, **backward**, **CrankNicolson**). The time step formula for **backward** time-scheme is for example :

.. code::

  deltaT = Foam::pow(3 x truncationError x Cmax[speciesi]/dC3dT3max[speciesi],1./3.)

where **dC3dT3maxmaximal** is the maximal value of the thid time derivative and **Cmax** the value of C in this cell.
