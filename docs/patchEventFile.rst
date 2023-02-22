.. _patchEventFile:

patchEventFile class
====================

Description
-----------

Event file handler/reader for the porousMultiphaseFoam toolbox which contains *m* dates with *n* flow rate values expressed in m3/s and attributed to a given patch name. Flow rate source terms are equally distributed over the surface of the given patch. Intermediate values are computed by linear interpolation between times. Flow rates before *time1* and after *timem* is equal to 0. 

Format
------

The class reads the file :

.. code::

    date time1
    patchName1 sourceTerm11
    patchName2 sourceTerm12
    ...
    patchNamen sourceTerm1n
    date time2
    patchName1 sourceTerm21
    patchName2 sourceTerm22
    ...
    date timem
    ...
    patchNamen sourceTermmn

Usage
-----

This class is used by :ref:`eventFlux` and :ref:`eventInfiltration`.
