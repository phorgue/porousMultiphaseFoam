.. _EulerD3dt3Scheme:

EulerD3dt3Scheme class
======================

Third-order Euler implicit d3dt3 using the current and three previous time-step values.

This class in mainly used by the :ref:`timestepManager` to evaluate truncation error error and compute simulation time step when second-order time scheme is used (*backward* for example).
