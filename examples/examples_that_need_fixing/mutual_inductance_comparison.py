#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:41:34 2019

@author: makinea1
"""

import numpy as np
import trimesh
from timeit import timeit

import pkg_resources
import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0,path)

from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix
from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix_from_A


#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/10x10_plane.obj')
coilmesh = trimesh.load(file_obj, process=False)

M1 = mutual_inductance_matrix(coilmesh, coilmesh)

M2 = mutual_inductance_matrix_from_A(coilmesh, coilmesh)

setup = (
        """
import numpy as np
import trimesh
from timeit import timeit

import pkg_resources
import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0,path)

from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix
from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix_from_A


#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/10x10_plane.obj')
coilmesh = trimesh.load(file_obj, process=False)
    """)

t1= timeit('M1 = mutual_inductance_matrix(coilmesh, coilmesh)', setup=setup, number=10)
t2 = timeit('M2 = mutual_inductance_matrix_from_A(coilmesh, coilmesh)', setup=setup, number=10)
