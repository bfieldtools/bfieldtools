# bfieldtools
A set of tools for magnetic field modelling with emphasis on magnetic field generated by continuous surface current density. The surface current density is decribed by a piecewise linear stream function on a triangle mesh.

The toolbox currently includes tools for
* Mesh operators (matrices) for mutual inductance, laplacian/resistance and gradient
* Magnetic field of current density on a mesh.
* Coil design 
* Thermal noise simulation 
* Spherical harmonics representation of data.
* Magnetic field of wire segments


# Dependencies
* [Numpy](https://www.numpy.org/)
* [Scipy](https://www.scipy.org/)
* [Trimesh](https://github.com/mikedh/trimesh)
* [CVXOPT](https://cvxopt.org/)
* [Quadpy](https://github.com/nschloe/quadpy/tree/master/quadpy)
* [Numba](https://numba.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [Mayavi](https://docs.enthought.com/mayavi/mayavi/)

# Documentation

For now, documentation is limited to function docstrings and comments. Some examples showcasing functionality are also included in the examples-folder.


## Installation instructions

Basic setuptools workflow applies.

`python setup.py build`
`python setup.py install`

*When installing an evolving piece of software as a package, take care where you install, which version etc.*


# License 

TBD.

