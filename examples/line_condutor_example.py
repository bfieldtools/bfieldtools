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
c.plot_loops()

#%% Calculate flux between two sets of line conductors
from bfieldtools.line_magnetics import flux

mesh2 = mesh.copy()
mesh2.vertices[:, 1] += 1
c2 = LineConductor(mesh=mesh2, scalars=sf)
f = flux(c, c2)

# Plot the mutual inductance matrix
import matplotlib.pyplot as plt

plt.matshow(f)
