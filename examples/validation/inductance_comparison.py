"""
Test diagonal of inductance matrix
===================================================

Use different number of quadrature points and two different implementations
"""

import numpy as np
import matplotlib.pyplot as plt

from bfieldtools.utils import load_example_mesh

coilmesh = load_example_mesh("10x10_plane")

from bfieldtools.mesh_impedance import self_inductance_matrix, mutual_inductance_matrix


M5 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=5)
M6 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=6)
M7 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=7)

plt.figure()
for m in (M5, M6, M7):
    plt.plot(np.diag(m))

MM5 = self_inductance_matrix(coilmesh, quad_degree=5, analytic_self_coupling=True)
MM6 = self_inductance_matrix(coilmesh, quad_degree=6, analytic_self_coupling=True)
MM7 = self_inductance_matrix(coilmesh, quad_degree=7, analytic_self_coupling=True)

#%% Plot the diagonals

plt.gca().set_prop_cycle(None)
for m in (MM5, MM6, MM7):
    plt.plot(np.diag(m), "--")
