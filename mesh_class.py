#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:41:43 2019

@author: Rasmus Zetter
"""

from mayavi import mlab
import trimesh


import utils
from laplacian_mesh import laplacian_matrix
from mutual_inductance_mesh import mutual_inductance_matrix


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

            self.verts = self.mesh.vertices
            self.tris = self.mesh.faces

            self.tri_areas = self.mesh.area_faces
            self.tri_normals = self.mesh.face_normals

        else: #Fall back on verts, tris parameters
            if isinstance(verts, type(None)) or isinstance(tris, type(None)):
                ValueError('You must provide either verts and tris or a mesh file')

            self.verts = verts
            self.tris = tris

            #Compute normals and surface already in constructor, should be quick
            self.tri_normals, self.tri_areas = utils.tri_normals_and_areas(self.verts, self.tris)



        #Useful mesh metrics etc, more can be added
        self.dual_areas = utils.dual_areas(self.tris, self.tri_areas)

        self.boundary_verts, self.inner_verts,
        self.boundary_tris, self.inner_tris = utils.find_mesh_boundaries(self.verts,
                                                                         self.tris,
                                                                         self.mesh.edges)
#        self.tri_barycenters = ...


        #Define uncomputed measures as None in constructor
        self.inductance = None
        self.resistance = None
        self.laplacian = None
        self.mass_matrix = None



    def compute_laplacian(self):
        '''
        Compute and return surface laplacian matrix as well as mass matrix.
        '''
        self.laplacian, self.mass_matrix = laplacian_matrix(self.verts, self.tris,
                                                            self.tri_normals,
                                                            self.tri_areas,
                                                            self.dual_areas)

        return self.laplacian, self.mass_matrix


    def compute_mutual_inductance(self):
        '''
        Compute and return mutual inductance matrix.
        '''

        self.inductance = mutual_inductance_matrix(self.verts, self.tris)

        return self.inductance


    def compute_resistance(self):
        '''
        Compute and return resistance/resistivity matrix.
        '''

        self.resistance = NotImplementedError('Function not implemented yet')

        return self.resistance

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



        return NotImplementedError()


if __name__ == '__main__':
    from matplotlib.tri import Triangulation
    import numpy as np

    xx = np.linspace(0, 1, 50)
    X, Y = np.meshgrid(xx, xx, indexing='ij')
    x = X.ravel()
    y = Y.ravel()
    z = np.zeros_like(x)
    print('Triangulating mesh')
    tt = Triangulation(x, y)

    verts = np.array([x, y, z]).T
    tris = tt.triangles


    obj = ToBeNamed(mesh_file='/m/nbe/project/hrmeg/matlab_koodit/CoilDesignPackage/CoilDesign/streamBEM/data/RZ_test_4planes_round.obj')

