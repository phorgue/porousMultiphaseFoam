.. _darcyGradPressureAniso:

darcyGradPressureAniso BC
=========================

Description
-----------
This boundary condition provides a fixed gra^dient condition for pressure **p** for a fixed velocity boundary condition when using the :ref:`anisoImpesFoam`. For more details about the pressure gradient computation, you can refer to :ref:`darcyGradPressure`.

In :ref:`anisoImpesFoam`, the permeability field **K** is *tensorial* and requires specific computation.

Usage
-----

.. code::

    <patchName>
    {
        type            darcyGradPressureAniso;
        phiName         phi; // optional (phi by default)
        value           uniform 0; // optional initial value (necessary for paraview display)
    }

Comments
--------

This boundary generally requires a *fixedValue* BC on velocity which may be variable in time and in space (**darcyGradPressureAniso** update pressure gradient at each time iteration).


Improvement
-----------

It is planned to merge this boundary condition to :ref:`darcyGradPressure` by developing a *templated* boundary conditions, accepting both scalar and tensorial permeability fields.
