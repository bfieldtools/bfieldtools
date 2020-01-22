#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:23:02 2020

@author: makinea1
"""


import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0, path)

import numpy as np
from bfieldtools.suhtools import SuhBasis
from trimesh.creation import icosphere
from mayavi import mlab
import trimesh
import matplotlib.pyplot as plt
from time import clock
from bfieldtools.utils import find_mesh_boundaries


# Create basis for a sphere (basis.eigenvals shere the same structure
# as spherical harmonic eigenvalues)
mesh = icosphere(4)

#import pkg_resources
##Load simple plane mesh that is centered on the origin
#file_obj = pkg_resources.resource_filename('bfieldtools',
#                'example_meshes/10x10_plane_hires.obj')
#mesh = trimesh.load(file_obj, process=True)
#t = np.eye(4)
#t[1:3,1:3] = np.array([[0,1],[-1,0]])
#mesh.apply_transform(t)

boundary, inner_verts = find_mesh_boundaries(mesh)

closed=True
basis = SuhBasis(mesh, 1, inner_vertices=inner_verts, closed_mesh=closed)

basis.Nc = 200
t0=clock()
basis.calculate_basis(closed_mesh=closed, shiftinvert=True)
print('Time with shift invert:', clock()-t0)
mlab.figure()
basis.plot(6)
e1= basis.eigenvals
b1=basis.basis
t0=clock()
basis.calculate_basis(closed_mesh=closed, shiftinvert=False)
print('Time without shift invert:', clock()-t0)
b2=basis.basis
mlab.figure()
basis.plot(6)
e2 = basis.eigenvals

plt.plot(e1)
plt.plot(e2)
