.. _eventFlux:

eventFlux BC
============

Description
-----------
This boundary condition provides a fixed value condition for a given variable **C** calculated as:

    **C = (constantValue+eventValue) / totalFluxBC**

where **constantValue** is a constant and **eventValue** a time-dependent fixed flux (linear interpolation between the two closest dates).

As for all source terms of the toolbox, negative values indicates that flow is going inside the domain.

Usage
-----

.. code::

    <patchName>
    {
        type            eventFlux;
        constantValue   -50;
        eventFile       "/path/to/event/file.dat";
        phiName         phi; // optional (phi by default)
        value           uniform 0; // optional initial value (necessary for paraview display)
    }

The event values are given in the event file **eventFile**.

Comments
--------

The keyword **eventFlux** is recognized by the dispersion models of the toolbox (such as **alphaDispersion**) which set the dispersion coefficient **Deff** equal to zero along the boundary **eventFlux** to ensure that fixed fluxes are correctly injected.

Restrictions
------------

The **eventFlux** BC can be used with:

- :ref:`groundwaterTransport2DFoam`
- :ref:`groundwaterTransportFoam`
- :ref:`porousScalarTransport2DFoam`
- :ref:`porousScalarTransportFoam`
