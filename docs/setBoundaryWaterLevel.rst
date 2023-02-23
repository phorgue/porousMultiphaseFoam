.. _setBoundaryWaterLevel:

setBoundaryWaterLevel utility
===============================

Description
-----------

Utility to set up head Pressure (3D) or potential (2D) on a specified patch

Can be used to set up water level for all groundwater solvers.

Usage
-----

1) to impose uniform head pressure (unsaturated solvers):

    .. code::

        setBoundaryWaterLevel -patch patchName -value 50.3

2) to impose uniform potential (saturated solvers):

    .. code::

        setBoundaryWaterLevel -patch patchName -value 50.3 -field potential

3) initialize pressure head with DEM  file (unsaturated solvers):

    .. code::

        setBoundaryWaterLevel -patch patchName -DEM stl_file

4) initialize potential with DEM  file (saturated solvers):

    .. code::

        setBoundaryWaterLevel -patch patchName -DEM stl_file -field potential

5)  initialize pressure head with STL  file (unsaturated solvers):

    .. code::

        setBoundaryWaterLevel -patch patchName -STL stl_file


Options
-------

`-patch` : **required** to specificy the patch to  set up water level

`-field` : **optional** <*h*> to specificy the field to initialize (*h* or *potential*)

`-value` : **optional** <*0*> for uniform initialization

`-DEM` : **optional** <*fileName*> file for DEM initialization

`-STL` : **optional** <*fileName*> file for STL initialization

`-threshold` : **optional** <*0*> minimum height for points to look in STL file

`-offset` : **optional** <*0*> specify the constant offset from the STL/DEM file

Restrictions
------------

- STL initialization can be only used for unsaturated solvers.
