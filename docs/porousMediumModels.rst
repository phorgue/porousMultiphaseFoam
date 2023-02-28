.. _porousMediumModels:

porousMediumModels
==================

This class allows the optional activation of a double porosity model for fluid flow in porous medium.

The porous medium model can be specified using the keyword :

.. code::

   porousMediumModel dualPorosity;

in the **constant/transportProperties** file.

If not defined, the solver considers the of the simple porosity model (*simplePorosity*).

.. warning::
   Dual porosity is available only for unsaturated Richards solver, i.e. :ref:`groundwaterFoam` and :ref:`groundwaterTransportFoam`.
   
Simple porosity
---------------

For the simple porosity model, you may specifiy uniform porosity or permeability in **constant/transportProperties**:

.. code::

   eps 0.305;
   K 1.25e-12;

Note that you may specify non-uniform permeability or porosity field in **constant/K** or **constant/eps**.

Dual porosity
-------------

In the case of dual porosity model, you may specify the properties of the second porous region in a dedicated dictionary *dualPorosityCoeffs* in the **constant/transportProperties** :

.. code::

   KFracture 1.1798e-12;

   dualPorosityCoeffs
   {
     a 0.01;
     beta 3;
     gammaW 0.4;
     Kmatrix 1.1798e-14;
     Kexchange 1.17982e-16;
   }


where :

- **a/beta/gammaW :**: the coefficient related to the dual porosity model
- **Kmatrix ::** the permeability of the matrix medium
- **Kexchange ::**  the permeability exchange between the two medium

Note that the *dualPorosity* model use **KFracture** field instead of the classical permeability **K** to avoid confusion  
  
As for simple porosity model, the permeabilities can be specified as non-uniform in **constant/KFracture**, **constant/Kmatrix** and **constant/Kexchange**.

Parameters for relative permeabilities and capillary pressure should be specified for each medium, *Fracture* and *Matrix*, in **constant/transportProperties**:

.. code::

   VanGenuchtenCoeffs
   {
     //- fracture properties
     thetaFracturemin 0;
     thetaFracturemax 0.025;
     alphaFracture 5.0; 
     mFracture 0.5;
     //- matrix properties
     thetaMatrixmin 0.10;
     thetaMatrixmax 0.475;
     alphaMatrix 0.5;
     mMatrix 0.333333;
   }
