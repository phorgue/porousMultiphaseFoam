.. _porousMediumTransportModels:

porousMediumTransportModels
===========================

This class handles the solving of scalar transport equations in porous media (advection/dispersion) for a simple or dual porosity model.

More information about the parametrization of the transport can be found in :ref:`porousScalarTransportFoam`.

Dual porosity
-------------

The dual porosity model for transport is linked to the dual porosity model for fluid flow (see :ref:`porousMediumModels`) which can be activated with:

.. code::

   porousMediumModel dualPorosity;

in the **constant/transportProperties** file.

If not defined, the solver considers the use of the simple porosity model (*simplePorosity*).

The transport dual porosity model uses the same parameters as for the flow (see :ref:`porousMediumModels`).


.. warning::
   Dual porosity model for transport is available only for unsaturated Richards solver, i.e. :ref:`groundwaterTransportFoam`.


