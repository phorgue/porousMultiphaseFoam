.. _index:

porousMultiphaseFoam Documentation
==================================

Welcome to the porousMultiphaseFoam **PMFv2406** documentation.

**PMF** is an open-source toolbox dedicated to simulation of flow and transport
processes in porous media. It is based on the OpenFOAM environment and therefore
benefits from its multiple feature such as pre- and post-processing tools or
parallel efficiency for example.

Initially developed for the modeling of two-phase flow in porous media (using
the IMPES method), PMF was then extended to the hydrology, with two
complementary approaches: a 3D unsaturated code solving the Richards' equation
and a low-cost 2D saturated code (using Dupuit's assumption).. 

A multi-specie transport module, coupled or decoupled, has been developed and
can be used either with unsaturated or saturated solvers. 

Multiple features have been added to help set up of realistic models :

- event file structures to handle time-varying source points or boundary
  conditions

- DEM/XY file structures to use topological data

- numerical methods for time-stepping and non-linear solvers.

Validation and development details can be found in this document or in the
:ref:`publications`


Content
=======

.. toctree::
   :maxdepth: 2

   installation
   solvers
   utilities
   libraries
   publications

.. toctree::
   :caption: Links

   Source repository <https://github.com/phorgue/porousMultiphaseFoam>
