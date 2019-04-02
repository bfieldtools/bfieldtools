#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 09:37:22 2019

@author: Antti Mäkinen
"""

import numpy as np
from scipy.sparse import csr_matrix, spdiags
from utils import tri_normals_and_areas, dual_areas


def laplacian_matrix(verts, tris):
    """
    Sparse Laplace(-Beltrami) operator

    Parameters:

        verts: (n, 3) array (float)
        tris: (m, 3) array (int) - indices into the verts array


    Returns:

        Cotangent weights: w_ij = - 0.5* (cot(alpha) + cot(beta))

    """
    N = verts.shape[0]
    R = verts[tris]  # Nt x 3 (corners) x 3 (xyz)
    # Edges opposite to the vertex
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2) # Nt x 3 (edges) x 3 (x,y,z)
    # Areas
    n, a = tri_normals_and_areas(verts, tris)
    ii = []
    jj = []
    cot = []
    # Loop over edges in triangles
    for i in range(3):
        i1 = (i+1) % 3
        i2 = (i+2) % 3
        ii.append(tris[:, i1])
        jj.append(tris[:, i2])
        cot.append(-0.5*(edges[:, i1, :]*edges[:, i2, :]).sum(axis=-1)/(2*a))

    ii = np.ravel(ii)
    jj = np.ravel(jj)
    cot = np.ravel(cot)
    # Build sparse matrix
    L = csr_matrix((cot, (ii, jj)), shape=(N, N), dtype=float)
    # Sum contribution from both triangles (alpha and beta angles)
    # neighbouring the edge
    L = L + L.T
    L = L - spdiags(L.sum(axis=0), 0, N, N)

    da = dual_areas(tris, a)
    # Area / Mass matrix
    A =spdiags(da, 0, verts.shape[0], verts.shape[0]).tocsr()

    return L, A


if __name__ == '__main__':
    """ Example for calculating and plotting eigen modes of Laplacian
    """
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    from scipy.linalg import eigh
    from mayavi import mlab

    xx = np.linspace(0, 1, 50)
    X, Y = np.meshgrid(xx, xx, indexing='ij')
    x = X.ravel()
    y = Y.ravel()
    z = np.zeros_like(x)
    print('Triangulating mesh')
    tt = Triangulation(x, y)

    verts = np.array([x, y, z]).T
    tris = tt.triangles
    L, M = laplacian_matrix(verts, tris)
    u, v = eigh(-L.todense(), M.todense())

    plt.plot(u)
    mlab.triangular_mesh(*verts.T, tris, scalars=v[:, 7])
