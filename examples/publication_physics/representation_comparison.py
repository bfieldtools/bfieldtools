#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 16:28:49 2019

@author: makinea1
"""

import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path in sys.path:
    sys.path.insert(0, path)


import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import trimesh
import pkg_resources

from bfieldtools.mesh_magnetics import magnetic_field_coupling
from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic
from bfieldtools.mesh_magnetics import scalar_potential_coupling
from bfieldtools.sphtools import compute_sphcoeffs_mesh, sphbasis
from bfieldtools.suhtools import suhbasis


# https://graphics.stanford.edu/~mdfisher/Data/Meshes/bunny.obj
mesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools',
                                                             'example_meshes/bunny_repaired.obj'), process=True)
# Bunny not still watertight, this fixes it
trimesh.repair.fill_holes(mesh)

mesh.vertices -= mesh.vertices.mean(axis=0)

mesh_field = mesh.copy()
mesh_field.vertices += 0.005*mesh_field.vertex_normals
mesh_field = trimesh.smoothing.filter_laplacian(mesh_field, iterations=1)

bsph = sphbasis(10)

Ca, Cb = bsph.basis_fields(mesh_field.vertices, 4)

bsuh = suhbasis(mesh, 20)
Csuh = magnetic_field_coupling_analytic(mesh, mesh_field.vertices) @ bsuh.basis

def plot_basis_fields(C, comps):
    d = 0.2
    i = 0
    j = 0
    for n in comps:
        p = 1.05*mesh_field.vertices.copy()
        p2 = mesh_field.vertices.copy()
#        p[:,1] -= i*d
#        p2[:,1] -= i*d
        p[:,0] += i*d
        p2[:,0] += i*d
        m = np.max(np.linalg.norm(C[:,:,n], axis=0))
        vectors = mlab.quiver3d(*p.T, *C[:,:,n].T, mode='arrow', colormap='Greys',
                                vmin=0, vmax=m)
        vectors.glyph.mask_input_points = True
        vectors.glyph.mask_points.maximum_number_of_points = 2000
        vectors.glyph.mask_points.random_mode_type = 1
        vectors.glyph.glyph_source.glyph_position = 'center'
        vectors.glyph.glyph_source.glyph_source.shaft_radius = 0.08
        vectors.glyph.glyph_source.glyph_source.tip_radius = 0.12
#        m = np.max(abs((C[:,:,n].T*mesh_field.vertex_normals.T).sum(axis=0)))
#        s =mlab.triangular_mesh(*p.T, mesh_field.faces,
#                             scalars=(C[:,:,n].T*mesh_field.vertex_normals.T).sum(axis=0),
#                             colormap='seismic', vmin=-m, vmax=m, opacity=0.7)
#        s.actor.property.backface_culling = True
        m = np.max(abs((C[:,:,n].T*mesh_field.vertex_normals.T).sum(axis=0)))
        s= mlab.triangular_mesh(*p2.T, mesh.faces,
                             scalars=(C[:,:,n].T*mesh_field.vertex_normals.T).sum(axis=0),
                             colormap='seismic', vmin=-m, vmax=m)
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.module_manager.scalar_lut_manager.number_of_colors = 32
        i+=1
#        if i==5:
#            i=0
#            j+=1

comps = [0,4,10,15]
scene = mlab.figure(bgcolor=(1,1,1), size=(1200, 300))
plot_basis_fields(np.moveaxis(Ca,2,0), comps)
scene.scene.parallel_projection = True
scene.scene.z_plus_view()
scene.scene.camera.position = [0.311195220798254, 0.015505198389291763, 1.5436867477960008]
scene.scene.camera.focal_point = [0.311195220798254, 0.015505198389291763, -0.009933453053236008]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 1.0, 0.0]
scene.scene.camera.clipping_range = [1.4018107339185428, 1.7484644679977115]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()

scene =mlab.figure(bgcolor=(1,1,1), size=(1200, 300))
plot_basis_fields(Csuh, comps)
#scene.scene.z_plus_view()
scene.scene.camera.position = [0.3064269505720481, -0.07122115917327232, 1.5412569063126864]
scene.scene.camera.focal_point = [0.311195220798254, 0.015505198389291763, -0.009933453053236008]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0010428412874734605, 0.998439994080673, 0.05582553808280774]
scene.scene.camera.clipping_range = [1.3838753615937698, 1.7677425271918352]
scene.scene.camera.position = [0.3064269505720481, -0.07122115917327232, 1.5412569063126864]
scene.scene.parallel_projection = True
