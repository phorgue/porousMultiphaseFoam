.. _XYfile:

XYfile class
============

Description
-----------

XY file handler/reader with linear interpolation (using three closest points).

Accepts random point distribution in file but it is much slower than DEMfile interpolation (due to
the 3 closest points search).

Format
------

The class reads the file :

.. code::

    x1 y1 value1
    x2 y2 value2
    ...
    xN yN valueN

Usage
-----

The class is used in the :ref:`setFieldsFromXY` for variable initialization based on a data file.
