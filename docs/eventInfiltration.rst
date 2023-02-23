.. _eventInfiltration:

eventInfiltration BC
====================

Description
-----------
This boundary condition provides a fixed value condition for the velocity variable **U** calculated as:

    **U = (constantValue+eventValue) & nf**

where **constantValue** is a constant, **eventValue** a time-dependent infiltration (linear interpolation between the two closest dates) and **nf** the normal to the patch.

As for all source terms of the toolbox, negative values indicates that flow is going inside the domain.

Usage
-----

.. code::

    <patchName>
    {
        type            eventInfiltration;
        constantValue   -50;
        eventFile       "/path/to/event/file.dat";
        value           uniform 0; // optional initial value (necessary for paraview display)
    }

The event values are given in the event file **eventFile**.

Restrictions
------------

The **eventInfiltration** BC can be used with:

- :ref:`anisoImpesFoam`
- :ref:`groundwater2DFoam`
- :ref:`groundwaterFoam`
- :ref:`groundwaterTransport2DFoam`
- :ref:`groundwaterTransportFoam`
- :ref:`impesFoam`
