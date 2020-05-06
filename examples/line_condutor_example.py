# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:30:12 2020

@author: Antti
"""

import numpy as np
import trimesh
from bfieldtools.line_conductor import LineConductor
from bfieldtools.mesh_conductor import StreamFunction
from bfieldtools.suhtools import SuhBasis
from bfieldtools.utils import load_example_mesh

mesh = load_example_mesh("10x10_plane")
sb = SuhBasis(mesh, 10)
sf = StreamFunction(sb.basis[:, 1], sb.mesh_conductor)
c = LineConductor(mesh=mesh, scalars=sf)


# Plot loops for testing
fig = c.plot_loops()

#%% Calculate mutual_inductance between two sets of line conductors

mesh2 = mesh.copy()
mesh2.vertices[:, 1] += 1
c2 = LineConductor(mesh=mesh2, scalars=sf)
c2.plot_loops(figure=fig)

Mself = c.line_mutual_inductance(c, separate_loops=True, radius=1e-2)
M2 = c.line_mutual_inductance(c2, separate_loops=True)


# Plot the mutual inductance matrix
import matplotlib.pyplot as plt

ff, ax = plt.subplots(1, 2)
plt.sca(ax[0])
plt.matshow(Mself, fignum=0)
plt.sca(ax[1])
plt.matshow(M2, fignum=0)
