# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 21:01:36 2019

@author: Antti
"""

#import sys
#path = 'C:/Users/Antti/Documents/Python Scripts/bfieldtools'
#if path not in sys.path:
#    sys.path.insert(0,path)

import numpy as np
from mayavi import mlab
from scipy.linalg import eigh
import trimesh
from bfieldtools import utils
from bfieldtools.laplacian_mesh import laplacian_matrix, mass_matrix, gradient

import pkg_resources

mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane.obj'))

mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)

L = laplacian_matrix(mesh)
M = mass_matrix(mesh)

u, v = eigh(-L.todense()[inner_verts][:,inner_verts], M.todense()[inner_verts][:,inner_verts])

#%%
from scipy.sparse import eye as speye
from scipy.sparse.linalg import spsolve
scalars = np.zeros(mesh.vertices.shape[0])
scalars[inner_verts] = v[:,45]
# Solve one-step of diffusion equation
#scalars2 = M @ spsolve(M - 0.001*L, scalars)
#scalars=scalars2
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=scalars)

#%%
N = 10
contours=np.linspace(scalars.min(),scalars.max(),N)
edge_vals = scalars[mesh.edges_unique]
contour_polys = []
for c in contours:
    # Get edges that contain the contour values
    edge_inds = (edge_vals.min(axis=1) <= c)*(edge_vals.max(axis=1) >= c)
    c_edge_vals = edge_vals[edge_inds]
    # Solve weights (barycentric coordinates) for each edge
    w0 = (c  - c_edge_vals[:,1])/(c_edge_vals[:,0]-c_edge_vals[:,1])
    w1 = (-c + c_edge_vals[:,0])/(c_edge_vals[:,0]-c_edge_vals[:,1])
    # Calculate points linearly interpolated from vertices
    points = mesh.vertices[mesh.edges_unique[edge_inds]]
    points = points[:,0]*w0[:, None] + points[:,1]*w1[:, None]
    
    # Determine adjacency
    c_edges_in_faces = edge_inds[mesh.faces_unique_edges]
    c_faces = np.any(c_edges_in_faces,axis=1)
    edges_in_c_faces = mesh.faces_unique_edges[c_faces]
    if len(edges_in_c_faces)==0:
        print('No contours at f=', c)
        continue
    # Each element of this array corresponds to a face containing the contour
    # The two values in the element are the edges (indices) adjacent to the face
    c_edges_in_c_faces = np.array([a[b] for a,b in zip(edges_in_c_faces, c_edges_in_faces[c_faces])])
    # Contour edges (indices pointing to the unique edge list)
    c_edge_inds = list(np.flatnonzero(edge_inds))
    # Init loop variables
    key=c_edges_in_c_faces[0,0]
    val=c_edges_in_c_faces[0,1]
    ii=0
    sorted_inds=[]
    kmax=len(c_edge_inds)
    k=0
    # Loop over c_edges´by essentially solving a linked list from c_edges_in_c_faces
    while k<kmax:
        sorted_inds.append(c_edge_inds.index(val))
        c_edges_in_c_faces[ii] =-1
        print(key)
        ii,jj = np.nonzero(c_edges_in_c_faces==key)
        if len(ii)==0:
            # Next edge not found in the adjacency list, contour must be closed now
            # Sort points containing contours by adjacency of the edges
            # and append to contour_polys
            contour_polys.append(points[sorted_inds])
            # Break the loop if all edges have been visited
            if np.all(c_edges_in_c_faces==-1):
                break
            # Else find a starting point in another contour at the same level
            sorted_inds=[]
            ii = np.flatnonzero(c_edges_in_c_faces[:,0] >= 0)[0]
            jj = 0
        else:
            # Edge found
            ii=ii[0]
            jj=jj[0]
        # Update key and value
        val = c_edges_in_c_faces[ii,jj]
        key = c_edges_in_c_faces[ii,(jj+1)%2]
        k+=1
        if k==kmax:
            raise RuntimeWarning('Something wrong with the contours, number of max iterations exceeded')
        
1

#%%
from scipy.sparse import spdiags
def simplify_contour(c, min_edge=1e-3, angle_threshold=2e-2, smooth=True):
    # Remove small edges by threshold
    vals = [np.ones(c.shape[0]), -np.ones(c.shape[0]), np.ones(c.shape[0])]
    D = spdiags(vals, [1,0,-c.shape[0]+1], c.shape[0], c.shape[0])
    edges =  D @ c
    c = c[np.linalg.norm(edges,axis=1) > min_edge]
    if len(c)==0:
        return None
    # Remove nodes on straigh lines
    D = spdiags(vals, [1,0,-c.shape[0]+1], c.shape[0], c.shape[0])
    H = spdiags(1/np.linalg.norm(D @ c, axis=1), 0, c.shape[0],c.shape[0])
    DD = H @ D
    c = c[np.linalg.norm(D.T @ DD @ c,axis=-1 ) > angle_threshold]
    if smooth:
            D = spdiags(vals, [1,0,-c.shape[0]+1], c.shape[0], c.shape[0])
            H = spdiags(1/np.linalg.norm(D @ c, axis=1), 0, c.shape[0],c.shape[0])
            DD = H @ D
            lengths = np.linalg.norm(D @ c, axis=1)
            lengths = 0.5*abs(D.T) @ lengths # Mean of edges
#            c = c - 0.2*lengths[:,None]*(D.T @ DD @ c)
            Nc = c.shape[0]
            c = spsolve(speye(Nc, Nc) + 1.0*spdiags(lengths,0,Nc,Nc)@(D.T @ DD), c)
    return c
    
contour_polys = [simplify_contour(c) for c in contour_polys]
contour_polys = [c for c in contour_polys if c is not None]

#%% Plot contours
for c in contour_polys:
    mlab.plot3d(*c[list(range(c.shape[0]))+[0]].T, color=(1,0,0), tube_radius=None)

#%%
#g = gradient(scalars, mesh, rotated=True)
#q = mlab.quiver3d(*mesh.triangles_center.T, *g, color=(0,0,1))
#q.glyph.glyph_source.glyph_position='center'
