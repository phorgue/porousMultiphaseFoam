.. _setFieldsFromXY:

setFieldsFromXY utility
=======================

Description
-----------

Initialize a given field y performing 3-points linear interpolation using
values given as table **x y value** (related to the class **XYfile**)

The file should be in the form :

.. code::

    x1 y1 value1
    x2 y2 value2
    ...
    xN yN valueN

with N the number of given points.

Usage
-----

.. code::

    setPermeabilityFieldFromXY -file inputFile_file -field fied_to_interpolate

Options
-------

`-file` : **required** to specificy the file used for 3-points linear interpolation

`-field` : **required** to specificy the field to interpolate

`-folder` : **optional** <*constant*> to specify the folder containing fields.

`-offset` : **optional** <*0*> to add an offset to the interpolated value

Restrictions
------------

- The interpolated field should be a **volScalarField**. 

- Only **x** and **y** values of the cell position are used for interpolation.

