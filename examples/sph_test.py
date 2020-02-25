#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:21:14 2020

@author: makinea1
"""


import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0, path)


import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

from bfieldtools.sphtools import basis_fields, ylm, derxlm, sinxlm

from trimesh.creation import icosphere
mesh = icosphere(1, 1)
print('Positions of evaluation points')
print(mesh.vertices)

theta=np.linspace(1e-6, np.pi*0.999)
sp = sinxlm(1,-1, theta)
plt.plot(sp)

aa, bb= basis_fields(mesh.vertices, 2)

for ii in range(aa.shape[1]):
    print('nans in sph component ', ii, ': ', np.isnan(aa[:,ii,:]).sum())

aa[np.isnan(aa)] = 100

# Plot the first component
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces,
                     scalars=np.linalg.norm(aa[:,0,:], axis=0))