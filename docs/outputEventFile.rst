.. _outputEventFile:

outputEventFile class
=====================

Description
-----------

Event file handler/reader for the porousMultiphaseFoam toolbox which contains *m* dates used for time output instead of classical openfoam write management.

Format
------

The class reads the file :

.. code::

    time1
    time2
    ...
    timen

Usage
-----

To use output event file, you should specify the file in **controlDict** as

.. code::

    outputEventFile "path/to/the/eventFile";

Restrictions
------------

The output file event can be used by all solvers.
