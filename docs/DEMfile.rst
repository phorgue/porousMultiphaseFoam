.. _DEMfile:

DEMfile class
=============

Description
-----------

DEM file handler/reader with bilinear interpolation fonction for the porousMultiphaseFoam toolbox

Format
------

The class reads the file :

.. code::

    x1 y1 z11
    x2 y1 z21
    ...
    xn y1 zn1
    x1 y2 z12
    x2 y2 z22
    ...
    xn y2 zn3
    x1 y3 z13
    ...
    xn ym znm

x and y values should be in ascending order with constant difference.

Usage
-----

In compatible solvers, the file name is specificied in **transportProperties** as

.. code::

  fileDEM "path/to/the/DEMFile";


Restrictions
------------

The DEM file can be used in:

- :ref:`groundwater2DFoam`
- :ref:`groundwaterTransport2DFoam`
- :ref:`setFieldsFromDEM`
- :ref:`subsetMeshFromDEM`
