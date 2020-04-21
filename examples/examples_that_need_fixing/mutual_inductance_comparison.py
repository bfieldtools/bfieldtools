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

from bfieldtools.utils import load_example_mesh

coilmesh = load_example_mesh("10x10_plane")

from bfieldtools.mesh_impedance import self_inductance_matrix, mutual_inductance_matrix

# M0 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=0)
# M1 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=1)
# M2 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=2)
# M3 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=3)
# M4 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=4)
M5 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=5)
M6 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=6)
M7 = mutual_inductance_matrix(coilmesh, coilmesh, quad_degree=7)

plt.figure()
for m in (M5, M6, M7):
    plt.plot(np.diag(m))


# MM1 = self_inductance_matrix(coilmesh, quad_degree=1)
# MM2 = self_inductance_matrix(coilmesh, quad_degree=2)
# MM3 = self_inductance_matrix(coilmesh, quad_degree=3)
# MM4 = self_inductance_matrix(coilmesh, quad_degree=4)
MM5 = self_inductance_matrix(coilmesh, quad_degree=5, analytic_self_coupling=True)
MM6 = self_inductance_matrix(coilmesh, quad_degree=6, analytic_self_coupling=True)
MM7 = self_inductance_matrix(coilmesh, quad_degree=7, analytic_self_coupling=True)


plt.gca().set_prop_cycle(None)
for m in (MM5, MM6, MM7):
    plt.plot(np.diag(m), "--")


from bfieldtools.mesh_impedance import (
    triangle_self_coupling,
    triangle_self_coupling_compact,
)

p1 = triangle_self_coupling(coilmesh)
p2 = triangle_self_coupling_compact(coilmesh)
