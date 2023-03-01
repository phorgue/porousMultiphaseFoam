.. _infiltrationEventFile:

infiltrationEventFile class
===========================

Description
-----------

Event file handler/reader for the porousMultiphaseFoam toolbox which contain *m* dates with one infiltration value expressed in m/s.

Intermediate values are computed by linear interpolation. Infiltration before *time1* and after *timem* are equal to 0.

.. note::
   You may specify an heterogenous time-variable infiltration by giving an infiltration value for each cell and each time. See tutorials groundwater2DFoam/timeVariableInfiltration for an example


Format
------

The class reads the file :

.. code::

    date time1
    infiltrationTerm1
    date time2
    infiltrationTerm2
    ...
    date timem
    infiltrationTermm

Usage
-----
To use source event file, you should specify the file in **transportProperties** as

.. code::

   infiltrationEventFile "flowRateOverTime.evt";


Restrictions
------------

The infiltration event file can be used in:

- :ref:`groundwater2DFoam`
- :ref:`groundwaterTransport2DFoam`
