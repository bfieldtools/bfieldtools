#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 12:16:53 2019

@author: Rasmus Zetter
"""

#%% Subdivide mesh, interpolate stream function to match finer mesh

from tvtk.api import tvtk

for stack in coil_axes:
    coil[stack].orig_I = coil[stack].I
    coil[stack].orig_mesh = coil[stack].mesh.copy()

    ug = tvtk.UnstructuredGrid(points=coil[stack].mesh.vertices)

    ug.set_cells(tvtk.Triangle().cell_type, coil[stack].mesh.faces)
    ug.point_data.scalars = coil[stack].I
    ug.point_data.scalars.name = 'scalars'


    coil[stack].mesh = coil[stack].orig_mesh.subdivide().subdivide()

    coil[stack].I =mlab.pipeline.probe_data(ug, *coil[stack].mesh.vertices.T)