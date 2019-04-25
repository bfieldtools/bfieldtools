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
    return M.T

@jit
def assemble_matrix2(tris1, tris2, Nverts1, Nverts2, triangle_data):
    """ Optimized  assembly of finite element matrix for
        precomputed triangle data for separate meshes 1 and 2

        Sums the triangle_data [Ntris (1), Ntris (2), 3 (nodes 1),3 (nodes 2)]
        for the nodes neighbouring the triangle
    """
    M = np.zeros((Nverts1, Nverts2))
    for i in range(tris2.shape[0]):  # Eval triangles
        for j in range(tris1.shape[0]):  # Source triangles
            for k in range(tris2.shape[1]): # Eval triangle hats
                for l in range(tris1.shape[1]): # Source triangle hats
                    M[tris2[i,k], tris1[j,l]] += triangle_data[i,j,k,l]
    return M.T

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



def fibonacci_sphere(samples=10, center=np.array([0, 0, 0]), radius = 1, randomize=True):
    '''
    Generates a set of points approximately evenly distributed on a sphere,
    with adjustable center and radius. Uses spherical Fibonacci Lattice.
    '''

    rnd = 1.
    if randomize:
        rnd = np.random.random() * samples

    points = np.zeros((samples, 3))
    offset = 2. / samples
    increment = np.pi * (3. - np.sqrt(5.));

    for i in range(samples):
        points[i, 1] = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - pow(points[i, 1], 2))

        phi = ((i + rnd) % samples) * increment

        points[i, 0] = np.cos(phi) * r
        points[i, 2] = np.sin(phi) * r

    return radius * points + center


def cylinder_points(radius=1, length=1, nlength=10, alpha=360, nalpha=10, center=np.array([0,0,0]), orientation=np.array([1,0,0])):
    '''
    Generate and return a set of points on a cylindrical surface.
    '''
    #Create the length array
    I = np.linspace(0, length, nlength)

    #Create alpha array avoid duplication of endpoints
    #Conditional should be changed to meet your requirements
    if int(alpha) == 360:
        A = np.linspace(0, alpha, num=nalpha, endpoint=False)/180*np.pi
    else:
        A = np.linspace(0, alpha, num=nalpha)/180*np.pi

    #Calculate X and Y
    X = radius * np.cos(A)
    Y = radius * np.sin(A)

    #Tile/repeat indices so all unique pairs are present
    pz = np.tile(I, nalpha)
    px = np.repeat(X, nlength)
    py = np.repeat(Y, nlength)

    points = np.vstack(( pz, px, py )).T

    #Shift to center
    shift = np.array(center) - np.mean(points, axis=0)
    points += shift

    #Orient tube to new vector

    #Grabbed from an old unutbu answer
    def rotation_matrix(axis,theta):
        a = np.cos(theta/2)
        b,c,d = -axis*np.sin(theta/2)
        return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                         [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                         [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

    ovec = orientation / np.linalg.norm(orientation)
    cylvec = np.array([1,0,0])

    if np.allclose(cylvec, ovec):
        return points

    #Get orthogonal axis and rotation
    oaxis = np.cross(ovec, cylvec)
    rot = np.arccos(np.dot(ovec, cylvec))

    R = rotation_matrix(oaxis, rot)
    return points.dot(R)
