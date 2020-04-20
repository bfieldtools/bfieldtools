# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 09:52:45 2020

@author: Antti

Study the eigenvalue spectrum of discretize self-inductance operator
of a spherical mesh. Compare the spectrum to analytical solution
"""

import numpy as np
from bfieldtools.mesh_properties import self_inductance_matrix, mutual_inductance_matrix
from bfieldtools.mesh_calculus import mass_matrix
from bfieldtools.utils import load_example_mesh
from scipy.linalg import eigh
import matplotlib.pyplot as plt

mesh = load_example_mesh("unit_sphere")

Nvals = 200
uu1 = []
uu2 = []
M = mass_matrix(mesh)

for q in range(4):
    L = self_inductance_matrix(mesh, analytic_self_coupling=True, quad_degree=q + 1)
    uu, vv = eigh(L, M.toarray(), eigvals=(0, Nvals - 1))
    uu1.append(uu)

    # L = self_inductance_matrix(mesh, analytic_self_coupling=False,
    # quad_degree=q+1)
    L = mutual_inductance_matrix(mesh, mesh, quad_degree=q + 1)
    uu, vv = eigh(L, M.toarray(), eigvals=(0, Nvals - 1))
    uu2.append(uu)


#%%
""" 
 Spherical harmonics are the eigenfunctions of self-inductance operator
 The correct eigenvalues derived using Taulu 2005 Eqs. (22, 23, A1, A5, A6)
 By considering the normal component of the magnetic field produced by 
 a single Y_lm. The result is
 e = mu_0*(l*(l+1)/(2*l+1))/R
"""
R = 1  # np.linalg.norm(mesh.vertices[mesh.faces].mean(axis=1),axis=-1).mean()
mu0 = (1e-7) * 4 * np.pi
ll = np.array([l for l in range(20) for m in range(-l, l + 1)])
ll = ll[:Nvals]
evals = ll * (ll + 1) / (2 * ll + 1)


#%%
# plt.plot(evals, 'k')
for u in uu1:
    uu_scaled = u / mu0 * R
    plt.plot(abs(evals[1:] - uu_scaled[1:]) / evals[1:], "-")
    # plt.plot(uu_scaled)

plt.gca().set_prop_cycle(None)
for u in uu2:
    uu_scaled = u / mu0 * R
    plt.plot(abs(evals[1:] - uu_scaled[1:]) / evals[1:], "--")
