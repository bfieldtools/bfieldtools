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
from mutual_inductance_mesh import mutual_inductance_matrix


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
            self.mesh = trimesh.load(mesh_file)


        else: #Fall back on verts, tris parameters
            if isinstance(verts, type(None)) or isinstance(tris, type(None)):
                ValueError('You must provide either verts and tris or a mesh file')
            self.mesh = trimesh.Trimesh(verts, tris)

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

        inductance = mutual_inductance_matrix(self.verts, self.tris)

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

    obj = ToBeNamed(mesh_file='/l/bfieldtools/example_meshes/macqsimal_testcoils_lowres.obj')

    obj.inductance;
    obj.laplacian;

    #for millimeters to meters
    obj.mesh.apply_scale(0.001)
    obj.verts = obj.mesh.vertices

    n_points = 250
    radius = 0.00075
    center = np.array([0, 0, 0])
    target_points = fibonacci_sphere(n_points, radius=radius, center=center)

    obj.C = compute_C(obj.mesh, target_points)

    target_field = 1e-8*np.ones((n_points, ))
    I = optimize_streamfunctions(obj, target_field,
                                 target_axis=2,
                                 target_error={'on_axis':0.05, 'off_axis':0.05},
                                 laplacian_smooth=0)

    limit = np.max(np.abs(I))

    mlab.figure()
    s=mlab.triangular_mesh(*obj.verts.T, obj.tris,scalars=I)
    s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
#    mlab.points3d(*target_points.T)

    B_target = np.vstack((obj.C[:, :, 0].dot(I), obj.C[:, :, 1].dot(I), obj.C[:, :, 2].dot(I))).T

    mlab.quiver3d(*target_points.T, *B_target.T)



    z = np.linspace(0, 0.03, 51)

    x = y = np.zeros_like(z)

    line_points = np.vstack((x, y, z)).T

    line_C = compute_C(obj.mesh, r=line_points)

    B_line = np.vstack((line_C[:, :, 0].dot(I), line_C[:, :, 1].dot(I), line_C[:, :, 2].dot(I))).T

    plt.plot(z*1e3, np.linalg.norm(B_line, axis=1)*1e9)
    plt.ylabel('Field amplitude (a.u., could be nT)')
    plt.xlabel('Position on z-axis [mm]')


#

    xx = np.linspace(0, 0.015, 31)
    Y, Z = np.meshgrid(xx, xx, indexing='ij')

    y = Y.ravel()
    z = Z.ravel()
    x = np.zeros_like(y)

    plane_points = np.vstack((x, y, z)).T

#    mlab.points3d(*plane_points.T)

    plane_C = compute_C(obj.mesh, plane_points)

    B_plane = np.vstack((plane_C[:, :, 0].dot(I), plane_C[:, :, 1].dot(I), plane_C[:, :, 2].dot(I))).T

#
##    mlab.points3d(*plane_points.T)
    mlab.quiver3d(*plane_points.T, *B_plane.T)

    mlab.contour3d(*plane_points.T, np.linalg.norm(B_plane, axis=1))