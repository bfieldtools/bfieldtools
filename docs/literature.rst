Literature
==========

The functionality of bfieldtools is based on a stream-function discretization of surface current onto a triangle mesh.

For a more theoretical and in-depth overview of the physics and computations used in bfieldtools, please see

.. bibliography:: references.bib
   :list: bullet
   :style: mystyle
   :filter: key == "makinen"

For a more software-centred  of bfieldtools and an overview of coil design, please see

.. bibliography:: references.bib
   :list: bullet
   :style: mystyle
   :filter: key == "zetter"

Selected papers
^^^^^^^^^^^^^^^^

A selection of useful papers can also be found below:

General introduction to stream functions

.. bibliography:: references.bib
   :list: bullet
   :style: mystyle
   :filter: key % "peeren" or key == "zeven"
   
Coil-design papers using the similar methods as bfieldtools

.. bibliography:: references.bib
   :list: bullet
   :style: mystyle
   :filter: key == "poole" or key == "pis" or key == "bring"
   
On the calculation of Laplace-Beltrami eigenfunctions (Surface Harmonics)

.. bibliography:: references.bib
   :list: bullet
   :style: mystyle
   :filter: key == "levy" or key == "reuter"
   
Thermal noise

.. bibliography:: references.bib
   :list: bullet
   :style: mystyle
   :filter: key == "roth" or key == "uhlemann"

   
Related software
================
 - libigl_ (C++ library for geometry processing with Python bindings)
 - gptoolbox_ (MATLAB toolbox for geometry processing) 
 - shtools_ (Fortran and Python tools for working with spherical harmonics)
 - `gBringout/CoilDesign`_   (MATLAB functions for coil design using stream functions and convex optimization)

.. _libigl: https://libigl.github.io/libigl-python-bindings/ 

.. _gptoolbox: https://github.com/alecjacobson/gptoolbox

.. _shtools: https://shtools.oca.eu/shtools/public/

.. _gBringout/CoilDesign: https://github.com/gBringout/CoilDesign
