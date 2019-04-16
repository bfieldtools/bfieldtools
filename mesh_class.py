#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:41:43 2019

@author: Rasmus Zetter
"""

from mayavi import mlab
import trimesh

import utils
from laplacian_mesh import laplacian_matrix, mass_matrix
from mutual_inductance_mesh import self_inductance_matrix, mutual_inductance_matrix


class LazyProperty():
    '''
    Implementation of lazily loading properties, see
    http://blog.pythonisito.com/2008/08/lazy-descriptors.html
    '''

    def __init__(self, func):
        self._func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__

    def __get__(self, obj, klass=None):
        if obj is None:
            return None
        result = obj.__dict__[self.__name__] = self._func(obj)
        return result


class ToBeNamed:
    '''
    Class that is used for surface mesh field calculations, e.g. coil design.
    Computation functions are typically external functions that are called
    using method wrappers.

    The mesh surface should typically consists of a single contiguous surface,
    although multi-surface meshes should also work. To combine multiple objects
    of this class, use the MultiSurface class (ToBeImplemented).
    '''

    def __init__(self, verts=None, tris=None, mesh_file=None):

        if mesh_file: #First, check if mesh_file passed
            self.mesh = trimesh.load(mesh_file, process=False)


        else: #Fall back on verts, tris parameters
            if isinstance(verts, type(None)) or isinstance(tris, type(None)):
                ValueError('You must provide either verts and tris or a mesh file')
            self.mesh = trimesh.Trimesh(verts, tris, process=False)

        self.verts = self.mesh.vertices
        self.tris = self.mesh.faces

        self.tri_areas = self.mesh.area_faces
        self.tri_normals = self.mesh.face_normals


        #Useful mesh metrics etc, more can be added
        self.dual_areas = utils.dual_areas(self.tris, self.tri_areas)

        self.boundary_verts, self.inner_verts,\
        self.boundary_tris, self.inner_tris = utils.find_mesh_boundaries(self.verts,
                                                                         self.tris,
                                                                         self.mesh.edges)


    @LazyProperty
    def laplacian(self):
        '''
        Compute and return surface laplacian matrix as well as mass matrix.
        '''
        laplacian = laplacian_matrix(self.verts, self.tris,
                                     self.tri_normals,
                                     self.tri_areas)

        return laplacian


    @LazyProperty
    def mass(self):
        '''
        Compute and return mesh mass matrix.
        '''
        mass = mass_matrix(self.verts, self.tris, self.tri_areas, self.dual_areas)

        return mass


    @LazyProperty
    def inductance(self):
        '''
        Compute and return mutual inductance matrix.
        '''

        #If mesh corresponds of many submeshes, compute these separately to save memory
        if self.mesh.body_count > 1:


            inductance = np.zeros((len(self.mesh.vertices), len(self.mesh.vertices)))

            #Split into separate sub-bodies
            submeshes = self.mesh.split(only_watertight=False)

            n_submeshes = len(submeshes)

            #Determine how submesh vertex indices correspond to full mesh vertex indices
            vertex_lookup = []

            for i in range(n_submeshes):
                vertex_lookup.append([])
                for vert in submeshes[i].vertices:
                    vertex_lookup[i].append(np.where((self.mesh.vertices == vert).all(axis=1) == True)[0][0])

            #Loop through block matrix components
            for i in range(n_submeshes):
                for j in range(n_submeshes):
                    if i==j:
                        sub_block = self_inductance_matrix(submeshes[i].vertices, submeshes[i].faces)
                    else:
                        sub_block = mutual_inductance_matrix(submeshes[i].vertices, submeshes[i].faces,
                                                 submeshes[j].vertices, submeshes[j].faces)

                    #Assign to full matrix
                    inductance[np.asarray(vertex_lookup[i])[:,None], vertex_lookup[j]] = sub_block


#            #Fill in lower triangle
#            i_lower = np.tril_indices(inductance.shape[0], -1)
#            inductance[i_lower] = inductance.T[i_lower]
        else:
            inductance = self_inductance_matrix(self.verts, self.tris)

        return inductance


    @LazyProperty
    def resistance(self, resistivity=1.68*1e-8, thickness=1e-4):
        '''
        Compute and return resistance/resistivity matrix using Laplace matrix.
        Default resistivity set to that of copper.
        NB! For now, the resistivity and thickness values are set in stone.
        To continue using @LazyProperty, they could be moved into class attributes.
        Alternatively, this LazyProperty could be turned into a method.
        '''

        R = resistivity / thickness
        resistance = R * self.laplacian.todense()

        #Set boundary vertices to zero
        resistance[obj.boundary_verts, :][:, obj.boundary_verts] = 0

        return resistance


    def plot_mesh(self):
        '''
        Simply plot the mesh surface.
        '''

        mesh = mlab.triangular_mesh(*self.verts.T, self.tris,
                                    representation='wireframe', color=(0, 0, 0))

        return mesh


    def plot_eigenmodes(self, n_modes=8):
        '''
        Plot eigenmodes of surface currents
        '''

        from scipy.linalg import eigh

        M = 0.5 * (self.inductance + self.inductance.T)
        Min = M[self.inner_verts[None, :], self.inner_verts[:, None]]
        print('Calculating modes')

        R = np.array(self.laplacian.todense())
        Rin = R[self.inner_verts[None, :], self.inner_verts[:, None]]
        w, v = eigh(-Rin, Min)

        mlab.figure()
        scalars = np.zeros((self.verts.shape[0], ))

        limit = np.max(abs(v[:, 0]))

        for ii in range(n_modes):

            n = int(np.sqrt(n_modes))

            i = ii % n
            j = int(ii / n)

            print(i, j)

            #Offset modes on XY-plane
            x = self.verts[:, 0] + i * 1.1
            y = self.verts[:, 1] + j * 1.1
            z = self.verts[:, 2]

            scalars[self.inner_verts] = v[:, ii]

            s = mlab.triangular_mesh(x, y, z, self.tris, scalars=scalars)

            s.module_manager.scalar_lut_manager.number_of_colors = 16
            s.module_manager.scalar_lut_manager.data_range = np.array([-limit, limit])

            s.actor.mapper.interpolate_scalars_before_mapping = True

        return s


if __name__ == '__main__':

    import numpy as np
    from utils import fibonacci_sphere
    from bringout_core import compute_C
    from coil_optimize import optimize_streamfunctions

    obj = ToBeNamed(mesh_file='/l/bfieldtools/example_meshes/macqsimal_testcoils_midres.obj')



    #for millimeters to meters
    obj.mesh.apply_scale(0.001)
    obj.verts = obj.mesh.vertices


    obj.inductance = self_inductance_matrix(obj.verts, obj.tris)
    obj.laplacian

    n_points = 200
    radius = 0.00075
    center = np.array([0, 0, 0])
    target_points = fibonacci_sphere(n_points, radius=radius, center=center)

    obj.C = compute_C(obj.mesh, target_points)

    n_stray_points = 60
    stray_radius = 0.02
    stray_points = fibonacci_sphere(n_stray_points, radius=stray_radius, center=center)


    obj.strayC = compute_C(obj.mesh, stray_points)

    target_field = 1e-10*np.ones((n_points, ))
    I, sol = optimize_streamfunctions(obj, target_field,
                                 target_axis=2,
                                 target_error={'on_axis':0.01, 'off_axis':0.01, 'stray':0.01},
                                 laplacian_smooth=0)

    limit = np.max(np.abs(I))




    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(480, 480))
    mlab.clf()

#    s=mlab.triangular_mesh(*obj.verts.T, obj.tris,scalars=I)

    surface = mlab.pipeline.triangular_mesh_source(*obj.verts.T, obj.tris,scalars=I)

    windings = mlab.pipeline.contour_surface(surface, contours=20)

#    s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
#    mlab.points3d(*target_points.T)

    B_target = np.vstack((obj.C[:, :, 0].dot(I), obj.C[:, :, 1].dot(I), obj.C[:, :, 2].dot(I))).T

    mlab.quiver3d(*target_points.T, *B_target.T)



    z = np.linspace(0, 0.03, 51)

    x = y = np.zeros_like(z)

    line_points = np.vstack((x, y, z)).T

    line_C = compute_C(obj.mesh, r=line_points)

    B_line = np.vstack((line_C[:, :, 0].dot(I), line_C[:, :, 1].dot(I), line_C[:, :, 2].dot(I))).T

    plt.plot(z*1e3, np.linalg.norm(B_line, axis=1)/np.mean(np.abs(target_field)))
    plt.ylabel('Field amplitude (target field units)')
    plt.xlabel('Position on z-axis [mm]')
    plt.grid(True)


#
#
#    xx = np.linspace(-0.015, 0.02, 15)
#    yy = np.linspace(0, 0.02, 15)
#    zz = np.linspace(0, 0.02, 15)
#    X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')
#
#    x = X.ravel()
#    y = Y.ravel()
#    z = Z.ravel()
#
#    grid_points = np.vstack((x, y, z)).T
#
#    mlab.points3d(*grid_points.T)
#
#    grid_C = compute_C(obj.mesh, grid_points)
#
#    B_grid = np.vstack((grid_C[:, :, 0].dot(I), grid_C[:, :, 1].dot(I), grid_C[:, :, 2].dot(I))).T
#    B_grid_matrix = B_grid.reshape((15, 15, 15, 3))
#
#    B_grid_matrix_norm = np.linalg.norm(B_grid_matrix, axis=-1)
#
#    field = mlab.pipeline.vector_field(X, Y, Z, B_grid_matrix[:,:,:,0], B_grid_matrix[:,:,:,1], B_grid_matrix[:,:,:,2],
#                                  scalars=B_grid_matrix_norm, name='B-field')
#
#    field2 = mlab.pipeline.vector_field(X, Y, -Z, B_grid_matrix[:,:,:,0], B_grid_matrix[:,:,:,1], B_grid_matrix[:,:,:,2],
#                                  scalars=B_grid_matrix_norm, name='B-field2')
#
#
#    vectors = mlab.pipeline.vectors(field,
#                          scale_factor=(X[1, 0, 0] - X[0, 0, 0]),
#                          )
#
#
#    vectors.glyph.mask_input_points = True
#    vectors.glyph.mask_points.on_ratio = 2
#
#    vcp = mlab.pipeline.vector_cut_plane(field)
#    vcp.glyph.glyph.scale_factor=10*(X[1, 0, 0] - X[0, 0, 0])
#    # For prettier picture:
#    #vcp.implicit_plane.widget.enabled = False
#
#    Bt=np.mean(np.linalg.norm(B_target, axis=1))
#
#    iso = mlab.pipeline.iso_surface(field,
#                                    contours=[0.005*Bt, 0.01*Bt, 0.05*Bt,0.1*Bt, 0.5*Bt,0.9*Bt, 1.1*Bt,2.5*Bt, 5*Bt],
#                                    opacity=0.6,
#                                    colormap='viridis')
#
#    iso2 = mlab.pipeline.iso_surface(field2,
#                                    contours=[0.005*Bt, 0.01*Bt, 0.05*Bt,0.1*Bt, 0.5*Bt,0.9*Bt, 1.1*Bt,2.5*Bt, 5*Bt],
#                                    opacity=0.6,
#                                    colormap='viridis')
#
#    # A trick to make transparency look better: cull the front face
#    iso.actor.property.frontface_culling = True
#
##
###    mlab.points3d(*plane_points.T)
##    mlab.quiver3d(*grid_points.T, *B_grid.T)
#
#    mlab.contour3d(X, Y, Z, B_grid_matrix/np.mean(np.abs(target_field)))