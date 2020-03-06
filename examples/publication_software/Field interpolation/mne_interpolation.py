#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 18:24:27 2020

@author: makinea1
"""

import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0, path)

import numpy as np
from bfieldtools.mesh_class import Conductor, StreamFunction
from bfieldtools.suhtools import SuhBasis
from trimesh.creation import icosphere
from mayavi import mlab
import trimesh
import matplotlib.pyplot as plt
from time import clock
from bfieldtools.utils import find_mesh_boundaries

import mne

from mne.datasets import sample
data_path = sample.data_path()
fname = data_path + '/MEG/sample/sample_audvis-ave.fif'
# Reading
condition = 'Left Auditory'
evoked = mne.read_evokeds(fname, condition=condition, baseline=(None, 0),
                      proj=True)
evoked.pick_types(meg='mag')
evoked.plot(exclude=[], time_unit='s')

i0, i1 = evoked.time_as_index(0.08)[0], evoked.time_as_index(0.09)[0]
field = evoked.data[:,i0:i1].mean(axis=1)

# Read BEM for surface geometry and transform to correct coordinate system
import os.path as op
subject = 'sample'
subjects_dir = op.join(data_path, 'subjects')
bem_fname = op.join(subjects_dir, subject, 'bem',
                    subject + '-5120-5120-5120-bem-sol.fif')
bem = mne.read_bem_solution(bem_fname)


# Head mesh 0
# Innerskull mesh 2
surf_index = 2

trans_fname = op.join(data_path, 'MEG', 'sample',
                      'sample_audvis_raw-trans.fif')
trans0 = mne.read_trans(trans_fname)
R = trans0['trans'][:3,:3]
t = trans0['trans'][:3,3]
# Surface from MRI to HEAD
rr = (bem['surfs'][surf_index]['rr'] - t) @ R
# Surface from HEAD to DEVICE
trans1 = evoked.info['dev_head_t']
R = trans1['trans'][:3,:3]
t = trans1['trans'][:3,3]
rr = (rr - t) @ R

mesh = trimesh.Trimesh(rr, bem['surfs'][surf_index]['tris'])
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)

surf_index = 0

R = trans0['trans'][:3,:3]
t = trans0['trans'][:3,3]
# Surface from MRI to HEAD
rr = (bem['surfs'][surf_index]['rr'] - t) @ R
# Surface from HEAD to DEVICE
R = trans1['trans'][:3,:3]
t = trans1['trans'][:3,3]
rr = (rr - t) @ R
head = trimesh.Trimesh(rr, bem['surfs'][surf_index]['tris'])

mlab.triangular_mesh(*head.vertices.T, head.faces, color=(0.5,0.5,0.5), opacity=0.5)

mesh=head



# Sensor locations and directions in DEVICE coordinate system
p = np.array([ch['loc'][:3] for ch in evoked.info['chs']
             if ch['ch_name'][-1]=='1' and ch['ch_name'][:3]=='MEG'])
n = np.array([ch['loc'][-3:] for ch in evoked.info['chs']
             if ch['ch_name'][-1]=='1' and ch['ch_name'][:3]=='MEG'])
# Plot sensor locations and directions
mlab.quiver3d(*p.T, *n.T, mode='arrow')

#%% Fit the surface current for the auditory evoked response
c = Conductor(mesh_obj=mesh, basis_name='suh', N_suh=150)
M = c.mass
#B_sensors = np.sum(c.B_coupling(p) * n[:,:,None], axis=1)
B_sensors = np.einsum('ijk,ij->ik',c.B_coupling(p), n)
#a = np.linalg.pinv(B_sensors, rcond=1e-15) @ field
ss = np.linalg.svd(B_sensors @ B_sensors.T, False, False)

#reg_exps = [0.5, 1, 2, 3, 4, 5, 6, 7, 8]
reg_exps = [1]
plot_this = True
rel_errors = []
for reg_exp in reg_exps:
    _lambda = np.max(ss)*(10**(-reg_exp))
    # Laplacian in the suh basis is diagonal
    BB = B_sensors.T @ B_sensors + _lambda*(-c.laplacian)/np.max(abs(c.laplacian))
    a = np.linalg.solve(BB, B_sensors.T@field)
    #a = B_sensors.T @ np.linalg.solve(BB, field)
    s = StreamFunction(a, c)
    b_filt = B_sensors @ s

    rel_error = np.linalg.norm(b_filt - field)/np.linalg.norm(field)
    print('Relative error:', rel_error*100, '%')
    rel_errors.append(rel_error)

    if plot_this:
        mlab.figure()
        surf = s.plot(False)
        surf.actor.mapper.interpolate_scalars_before_mapping = True
        surf.module_manager.scalar_lut_manager.number_of_colors = 16

#plt.plot(reg_exps, rel_errors, '.-')

# Additional filtering
#sbasis = SuhBasis(mesh, 100)
#A = sbasis.basis
#a_filt = A@(A.T @ M @ a)
#s = StreamFunction(a_filt, c)
#mlab.figure()
#s.plot(True)

#%% Interpolate to the sensor surface
import pkg_resources
#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                'example_meshes/meg_helmet.obj')
helmet = trimesh.load(file_obj, process=True)
# Bring the surface roughly to the correct place
helmet.vertices[:,2] -= 0.05

# Reset coupling by hand
c.B_coupling.reset()
mlab.figure()
B_surf = np.sum(c.B_coupling(helmet.vertices) * helmet.vertex_normals[:,:,None], axis=1)
#vecs = c.B_coupling(helmet.vertices)
mlab.quiver3d(*p.T, *n.T, mode='arrow')
scalars = B_surf @ s
surf = mlab.triangular_mesh(*helmet.vertices.T, helmet.faces, scalars=scalars,
                     colormap='seismic')
surf.actor.mapper.interpolate_scalars_before_mapping = True
surf.module_manager.scalar_lut_manager.number_of_colors = 15
surf2 = s.plot(False)
surf2.actor.mapper.interpolate_scalars_before_mapping = True
surf2.module_manager.scalar_lut_manager.number_of_colors = 15

#mlab.figure()
#U_surf = c.U_coupling(helmet.vertices)
#scalars = U_surf @ s
#surf = mlab.triangular_mesh(*helmet.vertices.T, helmet.faces, scalars=scalars,
#                     colormap='seismic')
#surf.actor.mapper.interpolate_scalars_before_mapping = True
#surf.module_manager.scalar_lut_manager.number_of_colors = 15
#surf2 = s.plot(False)
#surf2.actor.mapper.interpolate_scalars_before_mapping = True
#surf2.module_manager.scalar_lut_manager.number_of_colors = 15


#%% Scalar  pontential on a plane
##Load simple plane mesh that is centered on the origin
#file_obj = pkg_resources.resource_filename('bfieldtools',
#                'example_meshes/10x10_plane_hires.obj')
#plane = trimesh.load(file_obj, process=True)
##t = np.eye(4)
##t[1:3,1:3] = np.array([[0,1],[-1,0]])
##mesh.apply_transform(t)
#plane.vertices *= 0.03
#
#scalars = c.U_coupling(plane.vertices).max(axis=1)
#vert_mask = abs(scalars) > np.max(abs(scalars)/10)
#face_index = np.nonzero(plane.faces_sparse.T @ vert_mask)[0]
#plane = plane.subdivide(face_index)
#
#scalars = c.U_coupling(plane.vertices).max(axis=1)
#vert_mask = abs(scalars) > np.max(abs(scalars)/5)
#face_index = np.nonzero(plane.faces_sparse.T @ vert_mask)[0]
#plane = plane.subdivide(face_index)
#
#scalars = c.U_coupling(plane.vertices) @ s
#vert_mask = abs(scalars) > np.max(abs(scalars)/3)
#face_index = np.nonzero(plane.faces_sparse.T @ vert_mask)[0]
#plane = plane.subdivide(face_index)
#
#scalars = c.U_coupling(plane.vertices) @ s
#inner = abs(c.U_coupling(plane.vertices).sum(axis=1)) >1e-15
#scalars[inner] *= -1
#m = np.max(abs(scalars))/1.5
#surf1 = mlab.triangular_mesh(*plane.vertices.T, plane.faces, scalars=scalars,
#                     colormap='bwr', vmin=-m, vmax=m)
#surf1.actor.mapper.interpolate_scalars_before_mapping = True
#surf1.module_manager.scalar_lut_manager.number_of_colors = 15
#surf2 = s.plot(False)
#surf2.actor.mapper.interpolate_scalars_before_mapping = True
#surf2.module_manager.scalar_lut_manager.number_of_colors = 15
##mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, color=(1,1,1))


#%% Calculate magnetic field in a box
Nvol=30
x = np.linspace(-0.125, 0.125, Nvol)
vol_points = np.array(np.meshgrid(x, x, x, indexing='ij')).reshape(3,-1).T
#mlab.points3d(*vol_points.T)

c.B_coupling.reset()
Bvol_coupling = c.B_coupling(vol_points, Nchunks=100, analytic=True)
s = StreamFunction(a, c)
#s = StreamFunction(a, c)
Bvol = Bvol_coupling@s

#%% Plot with streamlines
from bfieldtools.mesh_calculus import gradient

#mlab.quiver3d(*vol_points.T, *Bvol.T)
mlab.figure(bgcolor=(1,1,1))
vecs = mlab.pipeline.vector_field(*vol_points.T.reshape(3,Nvol,Nvol,Nvol),
                                  *Bvol.T.reshape(3,Nvol,Nvol,Nvol))
vecnorm = mlab.pipeline.extract_vector_norm(vecs)

seed_points = mesh.vertices[mesh.faces].mean(axis=1) - 0.01*mesh.face_normals
#c1 = Conductor(mesh_obj=mesh, basis_name='vertex')
seed_vals = c.basis@c.inductance @ s
seed_vals_grad = np.linalg.norm(gradient(seed_vals, c.mesh), axis=0)
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=abs(seed_vals)**2, colormap='viridis')
seed_vals = abs(seed_vals[mesh.faces].mean(axis=1))**2
seed_vals[seed_vals_grad > seed_vals_grad.max()/1.8]=0
Npoints= 500
seed_inds = np.random.choice(np.arange(len(seed_vals)), Npoints, False, seed_vals/seed_vals.sum())
seed_points = seed_points[seed_inds]
#mlab.points3d(*seed_points.T, scale_factor=0.001)
#seed_vals /= seed_vals.max()
#rands = np.random.rand(len(seed_vals))
#seed_points = seed_points[seed_vals > rands]

streams=[]

for pi in seed_points:
    streamline = mlab.pipeline.streamline(vecnorm,
                                          integration_direction='both',
                                          colormap='BuGn',
                                          seed_visible=False,
                                          seedtype='point')
    streamline.seed.widget.position = pi
    streamline.stream_tracer.terminal_speed = 3e-13
    streamline.stream_tracer.maximum_propagation = 0.1
    streamline.actor.property.render_lines_as_tubes = True
    streamline.actor.property.line_width = 4.0
    streams.append(streamline)


# Magnetic flux
#s2 = StreamFunction(c.inductance @ s, c)
#mlab.figure()
#surf2 = s2.plot(False)
#surf2.actor.mapper.interpolate_scalars_before_mapping = True
#surf2.module_manager.scalar_lut_manager.number_of_colors = 15

#mlab.figure()
#surf2 = s.plot(False)
#surf2.actor.mapper.interpolate_scalars_before_mapping = True
#surf2.module_manager.scalar_lut_manager.number_of_colors = 16


# Custom colormap with alpha channel
streamine = streams[0]
lut = streamline.module_manager.scalar_lut_manager.lut.table.to_array()
lut[:, -1] = np.linspace(0, 255, 256)
streamline.module_manager.scalar_lut_manager.lut.table = lut
streamline.module_manager.scalar_lut_manager.data_range = np.array([1.0e-13,  1.0e-12])


##
for streamline in streams:
    streamline.stream_tracer.terminal_speed = 1e-13
    streamline.seed.widget.hot_spot_size = 0.1
    streamline.stream_tracer.initial_integration_step = 0.01
    streamline.stream_tracer.minimum_integration_step = 0.1

sensors = mlab.quiver3d(*p.T, *n.T, mode='cylinder')
sensors.glyph.glyph_source.glyph_source.height = 0.1
sensors.actor.property.color = (0.5, 0.5, 0.5)
sensors.actor.mapper.scalar_visibility = False
sensors.glyph.glyph_source.glyph_source.resolution = 32
sensors.glyph.glyph.scale_factor = 0.03

grad_s = gradient(c.basis@s, mesh, rotated=True)
q= mlab.quiver3d(*(mesh.vertices[mesh.faces].mean(axis=1).T), *grad_s,
              colormap='viridis', mode='arrow')
sensors.glyph.glyph_source.glyph_source.shaft_radius = 0.05
mlab.triangular_mesh(*head.vertices.T, head.faces, color=(0.8,0.8,0.8), opacity=1.0)

#
##    streamline.seed.widget.enabled = False
#    streamline.actor.property.line_width = 3.0


