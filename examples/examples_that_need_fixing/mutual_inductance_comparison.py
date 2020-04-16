#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:41:34 2019

@author: makinea1
"""

import numpy as np
import trimesh
from timeit import timeit
from time import clock

import pkg_resources
import sys
import matplotlib.pyplot as plt

# Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename(
    "bfieldtools", "example_meshes/10x10_plane.obj"
)
coilmesh = trimesh.load(file_obj, process=False)

from bfieldtools.mesh_properties import self_inductance_matrix, self_inductance_matrix2

# M0 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=0)
# M1 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=1)
# M2 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=2)
# M3 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=3)
# M4 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=4)
M5 = self_inductance_matrix(coilmesh, quad_degree=5)
M6 = self_inductance_matrix(coilmesh, quad_degree=6)
# M7 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=7)

for m in (M5, M6):
    plt.plot(np.diag(m))


# MM1 = self_inductance_matrix2(coilmesh, quad_degree=1)
MM2 = self_inductance_matrix2(coilmesh, quad_degree=2)
MM3 = self_inductance_matrix2(coilmesh, quad_degree=3)
MM4 = self_inductance_matrix2(coilmesh, quad_degree=4)
MM5 = self_inductance_matrix2(coilmesh, quad_degree=5)
MM6 = self_inductance_matrix2(coilmesh, quad_degree=6)

for m in (MM2, MM3, MM4, MM5, MM6):
    plt.plot(np.diag(m), "--")


from bfieldtools.mesh_properties import (
    triangle_self_coupling,
    triangle_self_coupling_compact,
)

p1 = triangle_self_coupling(coilmesh)
p2 = triangle_self_coupling_compact(coilmesh)
