darcyFoam utility
=================

Description
-----------

Utility to solve one phase flow in porous media. Soution overwrites current velocity/pressure field.

Usage
-----

.. code::

    darcyFoam -phase phaseName

Options
-------

`-phase`: **optional** <*a*> to specify the name of velocity field (**Ua** by default)

Configuration files
-------------------

**constant/transportProperties :**

.. code::

    phase.phaseName
    {
        mu              mu [ 1 -1 -1 0 0 0 0 ] 1e-3; // viscosity
        rho             rho [ 1 -3 0 0 0 0 0 ] 1000; // density
    }

Required fields
---------------

- **0/p :** pressure field

- **0/UphaseName :** velocity field (*phaseName* is given in the command line, **U** if not precised)

- **constant/g :** gravity field

- **constant/K :** permeability field
