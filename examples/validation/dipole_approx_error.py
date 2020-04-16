"""
Triangle dipole density approximation error
==================================================

Compare the exact solution for the potential of dipolar density with magnitude of a
linear shape function on a triangle with two approximations.

"""


from bfieldtools.integrals import (
    potential_vertex_dipoles,
    triangle_potential_dipole_linear,
    potential_dipoles,
)

from bfieldtools.mesh_calculus import mass_matrix
from bfieldtools.utils import load_example_mesh

import matplotlib.pyplot as plt


import numpy as np
from scipy.sparse import csc_matrix

mesh = load_example_mesh("unit_disc")


z = np.linspace(0, 10, 101)[1:]
x = np.zeros_like(z)
y = np.zeros_like(z)

points = np.array([x, y, z]).T

R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]
R_v = points[:, None, :] - mesh.vertices[None, :, :]

vertex_areas = mass_matrix(mesh, lumped=True).diagonal()

approx_pot_v = potential_vertex_dipoles(R_v, mesh.vertex_normals, vertex_areas)

approx_pot_f = potential_dipoles(R, mesh.face_normals, mesh.area_faces)
exact_pot_f = triangle_potential_dipole_linear(R, mesh.face_normals, mesh.area_faces)

# Map faces -> vertices
Nf = len(mesh.faces)
Nv = len(mesh.vertices)
M0 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 0])), (Nf, Nv))
M1 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 1])), (Nf, Nv))
M2 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 2])), (Nf, Nv))
exact_pot_v = (
    exact_pot_f[:, :, 0] @ M0 + exact_pot_f[:, :, 1] @ M1 + exact_pot_f[:, :, 2] @ M2
)
approx_pot_fv = (
    approx_pot_f[:, :, 0] @ M0 + approx_pot_f[:, :, 1] @ M1 + approx_pot_f[:, :, 2] @ M2
)

RE_fv = np.abs((approx_pot_fv - exact_pot_v) / exact_pot_v) * 100
RE_v = np.abs((approx_pot_v - exact_pot_v) / exact_pot_v) * 100


plt.semilogy(z, RE_fv.max(axis=1), label="potential_dipoles")
plt.semilogy(z, RE_v.max(axis=1), label="potential_vertex_dipoles")
plt.legend()

plt.ylabel("Maximum relative error (%)")
plt.xlabel("Distance")
