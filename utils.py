#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 17:59:56 2019

@author: Antti Mäkinen
"""

import numpy as np
import quadpy
from numba import jit

def tri_normals_and_areas(r, tri):
    """ Get triangle normals and areas from vertices (r) and
        triangle indices (tri)
    """
    n = np.cross(r[tri[:, 1], :]-r[tri[:, 0], :],
                 r[tri[:, 2], :]-r[tri[:, 0], :], axis=1)
    a = np.linalg.norm(n, axis=1)
    n = n/a[:, None]
    a /= 2
    return n, a


def get_quad_points(verts, tris, method='SevenPoint', index=None):
    """ Get quad points and weights from quadrature rules implemented in
        quadpy

        Returns:

        w: weights (Nquad)
        qp: quadrature points in each triangle (Ntris, Nquad)
    """
    methods = [k for k in quadpy.triangle.__dict__.keys() if k[0].isupper()]
    if method in methods:
        try:
            rule = quadpy.triangle.__dict__[method]()
        except(TypeError) as error:
            if index is not None:
                rule = quadpy.triangle.__dict__[method](index)
            else:
                print('The method requires index (check quadpy documentation)')
                raise error
    else:
        raise ValueError('method: '+method+' not in the following')
        print(methods)
    x = rule.points
    w = rule.weights

    qp = np.zeros((tris.shape[0], len(w), 3))
    for i, t in enumerate(tris):
        p0 = verts[t[0]]
        p1 = verts[t[1]]
        p2 = verts[t[2]]
        B = np.array([p1-p0, p2-p0])

        qp[i] = x @ B + p0

    return w, qp


@jit
def assemble_matrix(tris, Nverts, triangle_data):
    """ Optimized  assembly of finite element matrix for
        precomputed triangle data

        Sums the triangle_data [Ntris (1), Ntris (2), 3 (nodes 1),3 (nodes 2)]
        for the nodes neighbouring the triangle
    """
    M = np.zeros((Nverts, Nverts))
    for i in range(tris.shape[0]):  # Eval triangles
        for j in range(tris.shape[0]):  # Source triangles
            for k in range(tris.shape[1]): # Eval triangle hats
                for l in range(tris.shape[1]): # Source triangle hats
                    M[tris[i,k], tris[j,l]] += triangle_data[i,j,k,l]
    return M

@jit
def dual_areas(tris, ta):
    """ Calculate (dual) areas for each node in inds

        Dual area == area summed over the neighbouring triangles divide by 3
    """
    areas = np.zeros(np.max(tris)+1)
    for i in range(tris.shape[0]):
        for j in range(tris.shape[1]):
            areas[tris[i,j]] += ta[i]

    return areas/3


def find_mesh_boundaries(verts, tris, edges):
    '''
    Finds the open boundaries of a mesh by finding the edges that only
    belong to a single triangle. Returns an index array of inner vertices
    and triangles that do not touch the outer boundary.
    Takes edge parameter for convenience.
    '''

    unique, unique_idx, unique_count = np.unique(np.sort(edges, axis=-1), axis=0,
                                                 return_index=True,
                                                 return_counts=True)

    #If edge only used in one triangle, it is a boundary edge
    boundary_edges = edges[unique_idx[np.where(unique_count==1)]]

    #Create index arrays for boundary vertices
    boundary_verts = np.unique(boundary_edges.flatten())
    inner_verts = np.delete(np.arange(0, len(verts)), boundary_verts)

    #Find triangles using boundary vertices
    boundary_tris = np.array([], dtype=np.int)
    for vert in boundary_verts:
        boundary_tris = np.append(boundary_tris, np.where(np.any(tris == vert, axis=-1) is True)[0])

    #Create index arrays for boundary triangles
    boundary_tris = np.unique(boundary_tris)
    inner_tris = np.delete(np.arange(0, len(tris)), boundary_tris)

    return boundary_verts, inner_verts, boundary_tris, inner_tris

