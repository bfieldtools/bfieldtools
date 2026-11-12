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


M1 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=1)
M2 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=2)
M3 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=3)

plt.figure()
for m in (M1, M2, M3):
    plt.plot(np.diag(m))

MM1 = self_inductance_matrix(coilmesh, quad_degree=1, analytic_self_coupling=True)
MM2 = self_inductance_matrix(coilmesh, quad_degree=2, analytic_self_coupling=True)
MM3 = self_inductance_matrix(coilmesh, quad_degree=3, analytic_self_coupling=True)

#%% Plot the diagonals

plt.gca().set_prop_cycle(None)
for m in (MM1, MM2, MM3):
    plt.plot(np.diag(m), "--")
