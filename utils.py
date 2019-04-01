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
