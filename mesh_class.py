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
#        self.tri_barycenters =
#        self.

        #Define uncomputed measures as None in constructor
        self.inductance = None
        self.resistance = None
        self.laplacian = None
        self.mass_matrix = None


    def find_boundary(self):
        '''
        Finds the open boundaries of a mesh by finding the edges that only
        belong to a single triangle. Returns an index array of inner vertices
        and triangles that do not touch the outer boundary.
        '''
        unique, unique_idx, unique_count = np.unique(np.sort(self.mesh.edges,axis=-1), axis=0,
                                                     return_index=True,
                                                     return_counts=True)

        boundary_edges = self.mesh.edges[unique_idx[np.where(unique_count==1)]]

        self.boundary_verts = np.unique(boundary_edges.flatten())
        self.inner_verts = np.delete(np.arange(0, len(self.verts)), self.boundary_verts)

        self.boundary_tris = np.array([], dtype=np.int)
        for vert in self.boundary_verts:
            self.boundary_tris = np.append(self.boundary_tris, np.where(np.any(self.tris==vert, axis=-1)==True)[0])

        self.boundary_tris = np.unique(self.boundary_tris)

        self.inner_tris = np.delete(np.arange(0, len(self.tris)), self.boundary_tris)



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

