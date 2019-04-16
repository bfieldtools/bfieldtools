#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 09:01:32 2019

@author: makinea1
"""

import numpy as np
from time import clock
from utils import (get_quad_points, tri_normals_and_areas,
                   assemble_matrix, assemble_matrix2, dual_areas)


def gamma0(R, reg=1e-13, symmetrize=True):
    """ Integrals over the edges of a triangle called gamma_0

       (line charge potentials)

        Parameters:

            R : (N, 3, 3) array of points (Neval, Nverts, xyz)

        Returns:

            The analytic integrals for each vertex/edge

            array (Neval, Nverts)

        NOTE: MAY NOT BE VERY PRECISE FOR POINTS DIRECTLY AT TRIANGLE
        EDGES.
    """
    diffs = np.roll(R, 1, -2) - np.roll(R, 2, -2)
    dotprods1 = np.sum(np.roll(R, 1, -2)*diffs, axis=-1)
    dotprods2 = np.sum(np.roll(R, 2, -2)*diffs, axis=-1)
    dn = np.linalg.norm(diffs, axis=-1)
    del diffs
    n = np.linalg.norm(R, axis=-1)
    # Regularize s.t. neither the denominator or the numerator can be zero
    # Avoid numerical issues directly at the edge
    res = np.log((np.roll(n, 2, -1)*dn + dotprods2 + reg)
               / (np.roll(n, 1, -1)*dn + dotprods1 + reg))

    # Symmetrize the result since on the negative extension of the edge
    # there's division of two small values resulting numerical instabilities
    # (also incompatible with adding the reg value)
    if symmetrize:
        res2 = -np.log((np.roll(n, 2, -1)*dn - dotprods2 + reg)
                     / (np.roll(n, 1, -1)*dn - dotprods1 + reg))
        res = np.where(dotprods1+dotprods2 > 0, res, res2)
    res /= dn
    return res


def omega(R):
    """ Calculate the solid angle of a triangles using
        A. Van Oosterom and J. Strackee
        IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING,
        VOL. BME-30, NO. 2, 1983

        Parameters:

            R : (N, ..., 3, 3) array of points (Neval, ..., Nverts, xyz)

        Points correspond to relative coordinates (x,y,z) of
        N triangles/evaluation points for
        the 3 corners of the triangles/triangle

        Neval can be number of evaluation points for the same triangle
        or number of triangles for the same evaluation points

        Returns:

            Solid angles of triangle(s) at evaluation points
    """
    # Distances
    d = np.linalg.norm(R, axis=-1)
    # Scalar triple products
    stp = np.linalg.det(R)
    # Denominator
    denom = np.prod(d, axis=-1)
    for i in range(3):
        j = (i+1) % 3
        k = (i+2) % 3
        denom += np.sum(R[..., i, :]*R[..., j, :], axis=-1)*d[..., k]
    # Solid angles
    sa = 2*np.arctan2(stp, denom)
    return sa


def triangle_potential(R, tn, planar = False):
    """ 1/r potential of a uniform triangle

        Parameters:

            R : (N, (Ntri), 3, 3) array of displacement vectors
               (Neval, ...., Ntri_verts, xyz)

               These are displacement vectors

            tn : ((Ntri), 3) array of points (Ntri, dir)
    """
    if len(R.shape) > 3:
        tn_ax = tn[:, None, :]
    else:
        tn_ax = tn
    summands = gamma0(R)*np.sum(tn_ax*np.cross(np.roll(R, 1, -2),
                                np.roll(R, 2, -2), axis=-1), axis=-1)
    if not planar:
        csigned = np.sum(np.take(R, 0, -2)*tn, axis=-1)
        result = np.sum(summands, axis=-1) - csigned*omega(R)
    else:
        print('Assuming all the triangles are in the same plane!')
        result = np.sum(summands, axis=-1)
    return result




def mutual_inductance_matrix(verts1, tris1, verts2, tris2,
                             tri_normals1=None, tri_areas1=None,
                             tri_normals2=None, tri_areas2=None,
                             planar=False):
    """ Calculate a mutual inductance matrix for hat basis functions
        (stream functions) in the triangular mesh described by

        verts1: Nv x 3 array of mesh1 vertices (coordinates)
        tris1: Nt x 3 array of mesh1 triangles (indices to verts array)

        verts2: Nv x 3 array of mesh2 vertices (coordinates)
        tris2: Nt x 3 array of mesh2 triangles (indices to verts array)
    """
    R = verts1[tris1]  # Nt x 3 (corners) x 3 (xyz)
    # Calculate quadrature points
    weights, quadpoints = get_quad_points(verts2, tris2, 'Centroid')
    # Nt x Nquad x  3 (x,y,z)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_normals1, type(None)) or isinstance(tri_areas1, type(None)):
        tri_normals1, tri_areas1 = tri_normals_and_areas(verts1, tris1)

    if isinstance(tri_normals2, type(None)) or isinstance(tri_areas2, type(None)):
        tri_normals2, tri_areas2 = tri_normals_and_areas(verts2, tris2)

    RR = quadpoints[:,:, None, None, :] - R[None, None, :, :, :]
    print('Calculating potentials')
    pots = triangle_potential(RR, tri_normals1, planar=planar) # Ntri_eval, Nquad, Ntri_source
    pots = np.sum(pots*weights[None,:,None], axis=1) # Ntri_eval, Ntri_source

    # Calculate edge vectors for each triangle
    edges1 = np.roll(R, 1, -2) - np.roll(R, 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)
    edges2 = np.roll(verts2[tris2], 1, -2) - np.roll(verts2[tris2], 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)
    tri_data = np.sum(edges1[None,:,None,:,:]*edges2[:,None,:,None,:], axis=-1) # i,j,k,l
    tri_data /= (tri_areas2[:,None]*tri_areas1[None,:]*4)[:,:,None,None]
    tri_data *= (tri_areas2[:, None]*pots)[:,:,None,None]
    print('Inserting stuff into M-matrix')

    M = assemble_matrix2(tris1, tris2, verts1.shape[0], verts2.shape[0], tri_data)
    return M*1e-7

def self_inductance_matrix(verts, tris, tri_normals=None, tri_areas=None, planar=False):
    """ Calculate a self inductance matrix for hat basis functions
        (stream functions) in the triangular mesh described by

        verts: Nv x 3 array of mesh vertices (coordinates)
        tris: Nt x 3 array of mesh triangles (indices to verts array)
    """
    R = verts[tris]  # Nt x 3 (corners) x 3 (xyz)
    # Calculate edge vectors for each triangle
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)
    # Calculate quadrature points
    weights, quadpoints = get_quad_points(verts, tris, 'Centroid')
    # Nt x Nquad x  3 (x,y,z)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_normals, type(None)) or isinstance(tri_areas, type(None)):
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)

    RR = quadpoints[:,:, None, None, :] - R[None, None, :, :, :]
    print('Calculating potentials')
    pots = triangle_potential(RR, tri_normals, planar=planar) # Ntri_eval, Nquad, Ntri_source
    pots = np.sum(pots*weights[None,:,None], axis=1) # Ntri_eval, Ntri_source


    tri_data = np.sum(edges[None,:,None,:,:]*edges[:,None,:,None,:], axis=-1) # i,j,k,l
    tri_data /= (tri_areas[:,None]*tri_areas[None,:]*4)[:,:,None,None]
    tri_data *= (tri_areas[:, None]*pots)[:,:,None,None]
    print('Inserting stuff into M-matrix')

    # Non optimized version of the matrix assembly
#    M = np.zeros((verts.shape[0], verts.shape[0]))
#    for i, t1 in enumerate(tris):  # Eval triangles
#        for j, t2 in enumerate(tris):  # Source triangles
#            for k, v1 in enumerate(t1): # Eval triangle hats
#                for l, v2 in enumerate(t2): # Source triangle hats
##                    e1 = edges[i, k]
##                    e2 = edges[j, l]
##                    gradproduct = np.sum(e1*e2)/(a[i]*a[j]*4)
##                    M[v1, v2] += pots[i, j]*a[i]*gradproduct
    M = assemble_matrix(tris, verts.shape[0], tri_data)
    return M*1e-7




if __name__ == '__main__':
    """ Example for calculating independent current modes in a wall with
        a door
    """
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    from scipy.linalg import eigh
    from laplacian_mesh import laplacian_matrix

    xx = np.linspace(0, 1, 50)
    X,Y = np.meshgrid(xx,xx,indexing='ij')
    x = X.ravel()
    y = Y.ravel()
    z = np.zeros_like(x)
    print('Triangulating mesh')
    tt = Triangulation(x,y)

    verts = np.array([x,y,z]).T
    tris = tt.triangles

    # Determine inner triangles, this is a mess currently
    # Could be done in a systematic manner
    centers = np.mean([x[tris], y[tris]], axis=-1)
    trimask = ~((centers[1] < 0.7)*(abs(centers[0]-0.2) < 0.1))
    tris2 = tris[trimask]
    tt2 = Triangulation(x,y, tris2)
    boundary_tris = np.array([-1 in n for n in tt2.neighbors])
    all_inds = np.unique(tris.flatten())
    boundary_inds0 = np.unique(tris2[boundary_tris]).flatten()
    boundary_inds = [b for b in boundary_inds0 if np.sum(tris2==b)<=4]
    boundary_inds = np.array(list(boundary_inds) +
                          list(np.nonzero((abs(y-0.7)<0.01)*(abs(x-0.3)<0.01))[0]))
#    plt.triplot(x, y, tris)
#    plt.plot(x[boundary_inds], y[boundary_inds], 'r*')

    print('Calculating triangle stuff')
    n, a = tri_normals_and_areas(verts, tris)
    da = dual_areas(tris, a)

    boundary_all = (np.isclose(x,0))+(np.isclose(y,0))+(np.isclose(x,1))+(np.isclose(y,1)) > 0

    inner_inds = np.setdiff1d(all_inds,boundary_inds)
    inner_inds = np.setdiff1d(inner_inds,np.nonzero(boundary_all)[0])


    M = self_inductance_matrix(verts, tris)
    M=0.5*(M+M.T)
    Min = M[inner_inds[None,:], inner_inds[:,None]]
    print('Calculating modes')
    L = laplacian_matrix(verts, tris)
    L = np.array(L.todense())
    w,v = eigh(-L[inner_inds[None,:], inner_inds[:,None]], Min)

#%% Plot eigenmodes of surface currents on thin wall
    from mayavi import mlab
    mlab.figure()
    scalars = np.zeros(x.shape)
    Nmodes = 16
    limit = np.max(abs(v[:,0]))
    for ii in range(Nmodes):
        n = int(np.sqrt(Nmodes))
        i = ii % n
        j = int(ii/n)
        print(i,j)
        x = verts[:,0] + i*1.1
        y = verts[:,1] + j*1.1
        z = verts[:,2]
        scalars[inner_inds] = v[:,ii]
#        scalars[inner] = v[:,4] +v[:,5]
        s=mlab.triangular_mesh(x,y,z, tris, scalars=scalars) #M[:,70])
        s.module_manager.scalar_lut_manager.number_of_colors = 16
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True

