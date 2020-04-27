"""
Computation time for shiftinvert eigenvalue decomp
===========================================================

"""


import numpy as np
from bfieldtools.suhtools import SuhBasis
from mayavi import mlab
import trimesh
import matplotlib.pyplot as plt
from time import clock


# Create basis for a sphere (basis.eigenvals shere the same structure
# as spherical harmonic eigenvalues)
mesh = trimesh.creation.icosphere(4)

closed = False
basis = SuhBasis(mesh, 1)

# Choose Nc and recalculate basis with shift-invert and without
basis.Nc = 100
t0 = clock()
basis.calculate_basis(shiftinvert=True)
print("Time with shift invert:", clock() - t0)
f = mlab.figure()
basis.plot(15, figure=f)
e1 = basis.eigenvals
b1 = basis.basis
t0 = clock()
basis.calculate_basis(shiftinvert=False)
print("Time without shift invert:", clock() - t0)
b2 = basis.basis
f = mlab.figure()
basis.plot(15, figure=f)
e2 = basis.eigenvals

plt.plot(e1)
plt.plot(e2)
