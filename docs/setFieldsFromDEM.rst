.. _setFieldsFromDEM:

setFieldsFromDEM utility
========================

Description
-----------

Set the permeability field by performing bilinear interpolation using
DEM file which should be in the form :

.. code::

    x1 y1 z11
    x2 y1 z21
    ...
    xn y1 zn1
    x1 y2 z12
    x2 y2 z22
    ...
    xn zm znm

with x sorted in ascending order. delta_x and delta_y should be constant.

Usage
-----

.. code::

    setFieldsFromDEM -file DEMfile -field field_to_initialize

Options
-------

`-file` : **required** to specificy the DEM file used for bi-linear interpolation

`-field` : **required** to specificy the field to interpolate

`-folder` : **optional** <*constant*> to specify the folder containing field.

`-offset` : **optional** <*0*> to add an offset to the interpolated value

Restrictions
------------

- the interpolated field should be a **volScalarField**
