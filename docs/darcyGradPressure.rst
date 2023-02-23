.. _darcyGradPressure:

darcyGradPressure BC
====================

Description
-----------
This boundary condition provides a fixed gradient condition for pressure **p**, pressure head **h** or potential **potential** to impose fixed Darcy velocity. After computing the flux **phi** corresponding to the fixed velocity **U**, fixed gradient is computed as:

    **grad(p) = - (phi-phiG-phiPc)/M**

where **phi** is the total flux, **phiG** the gravity-induced flux, **phiPc** capillary-induced flux and **M** the phase mobility (depends on the solver used)

Solver specificities
--------------------

- all solvers except :ref:`impesFoam`: **phiPc = 0**

- saturated *quasi-2D* solvers (:ref:`groundwater2DFoam`, :ref:`groundwaterTransport2DFoam`) : **phiG = 0**

- Phase mobility expression:

    - Two-phase flow solver: **M = K x kr / mu** with **K** the permeability, **kr** the relative permeability and **mu** the dynamic viscosity.

        :ref:`impesFoam`

    - saturated *quasi-2D* solvers: **K x g x rho / mu**: with **g** the gravity and **rho** the density

        :ref:`groundwater2DFoam`, :ref:`groundwaterTransport2DFoam`

    - unsaturated *groundwater* solvers: **K x g x rho x kr / mu**

        :ref:`groundwaterFoam`, :ref:`groundwaterTransportFoam`


Usage
-----

.. code::

    <patchName>
    {
        type            darcyGradPressure;
        phiName         phi; // optional (phi by default)
        value           uniform 0; // optional initial value (necessary for paraview display)
    }

Comments
--------

This boundary generally requires a *fixedValue* BC on velocity which may be variable in time and in space (**darcyGradPressure** updates pressure gradient at each time iteration).

Restrictions
------------

This boundary works for most solvers except :ref:`anisoImpesFoam` which requires the use of  :ref:`darcyGradPressureAniso` due to its *tensorial* permeability field.
