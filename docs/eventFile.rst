.. _eventFile:

eventFile class
===============

Description
-----------

The is an general class used to handle/read the various event files for the porousMultiphaseFoam toolbox. This can be used to impose volumic or mass flow rate source inside the domain or on specified patch. It can be used to impose time variable uniform infiltration or specify user-defined time outputs.

Tips
----

If you want to impose step function with events, you may define two identical dates with different values. For example, if we want to get an infiltration equal to 11 from 10 seconds to 20 seconds and then equal to 22 from 20 seconds to 30, the eventInfiltration file gives:

.. code::

    date 10
    11
    date 20
    11
    date 20
    22
    date 30
    22

.. toctree::
   :maxdepth: 1

   outputEventFile
   patchEventFile
   sourceEventFile
   infiltrationEventFile
