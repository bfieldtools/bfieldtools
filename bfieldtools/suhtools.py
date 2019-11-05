#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:02:30 2019

@author: makinea1

Class for calculating the surface harmonic representation of
the magnetic field. Surface harmonics can represent
of any stream function as a series.

Surface harmonics == Laplace-Beltrami eigenfunctions

"""

from bfieldtools.laplacian_mesh import laplacian_matrix, mass_matrix
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.mesh_class import CouplingMatrix
from scipy.sparse.linalg import eigsh
import numpy as np

class suhbasis():
    """ Class for representing magnetic field using surface harmonics
    """

    def __init__(self, mesh, Nc, closed_mesh=False):
        """
        Parameters
            mesh : Trimesh-object representing the boundary on which
                        current density is specified

            Nc : Number of components
        """

        self.mesh = mesh
        self.Nc = Nc
        self.calculate_basis(closed_mesh)


    def calculate_basis(self, closed_mesh=True):
        """ Calculate basis functions as eigenfunctions of the laplacian
        
            closed_mesh: if True, calculate the basis for the whole mesh
                         if False, calculate for inner vertices (zero-dirichlet condition)
                                    and use zero for the boundary 
        """

        if closed_mesh:
            assert self.mesh.is_watertight

        L = laplacian_matrix(self.mesh)
        M = mass_matrix(self.mesh)

        if closed_mesh:
            v0 = np.ones(L.shape[1]) # avoid random basis for symmetric geometries
            vals, funcs = eigsh(-L, self.Nc+1, M, which='SA', v0 = v0)
    
            # The first function is constant and does not yield any field
            self.basis = funcs[:,1:]
            self.eigenvals = vals[1:]
        else:
            from .utils import find_mesh_boundaries
            self.boundary_verts, \
            self.inner_verts,\
            self.boundary_tris, \
            self.inner_tris = find_mesh_boundaries(self.mesh.vertices, 
                                                   self.mesh.faces, self.mesh.edges)
            Linner = L[self.inner_verts][:, self.inner_verts]
            Minner = M[self.inner_verts][:, self.inner_verts]
            v0 = np.ones(Linner.shape[1]) # avoid random basis for symmetric geometries
            vals, funcs = eigsh(-Linner, self.Nc, Minner, which='SA', v0 = v0)
    
            # Insert values to the inner verts
            self.basis = np.zeros((L.shape[0], self.Nc))
            self.basis[self.inner_verts] = funcs
            self.eigenvals = vals[1:]
            

    def field(self, coeffs, points):
        """ Calculate field at points

            Parameters

            coeffs : (Nc,) array of basis function coefficients
            points : (N_points, 3) field evaluation points

            Returns:

                field : (N_points, 3) magnetic field
        """
        return (self.basis_fields(points) @ coeffs).T

    def basis_fields(self, points):
        """ Calculate basis fields at points

            Return:

                Fields (3, N_points, N_coeffs)
        """
        C = compute_C(mesh, points)

        return C @ self.basis

    def fit_coeffs(self, points, data):
        """ Fit basis function coefficients to the data
        """
        assert len(data) > self.Nc

        A = self.basis_fields(points).reshape(-1, self.Nc)
        b = data.T.flatten()
        x, res, rank, s = np.linalg.lstsq(A, b, rcond=None)

        if rank < self.Nc:
            print('Matrix rank not full, result might be in accurate')
        return  x



if __name__ == '__main__':
    """ Simple testing script
    """

    from trimesh.creation import icosphere
    from mayavi import mlab

    # Create basis for a sphere (basis.eigenvals shere the same structure
    # as spherical harmonic eigenvalues)
    mesh = icosphere(4)
    basis = suhbasis(mesh, 20)

    s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces,
                             scalars=basis.basis[:,6])
    s.enable_contours = True


    coeffs = np.zeros(basis.basis.shape[1])
    coeffs[6] = 1

    # Plot field outside
    points = mesh.vertices * 2.0
    field = basis.field(coeffs, points)
    v = mlab.quiver3d(*points.T, *field.T)
    v.glyph.glyph_source.glyph_position = 'center'

    # Plot field inside
    points = mesh.vertices * 0.5
    field = basis.field(coeffs, points)
    v = mlab.quiver3d(*points.T, *field.T)
    v.glyph.glyph_source.glyph_position = 'center'

    # Modify the field and fit coeffs
    field1 = np.zeros_like(field)
    field1[:, 1] = np.max(field)
    field2 = field + field1

    coeffs2 = basis.fit_coeffs(points, field2)

    mlab.figure()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=basis.basis @ coeffs2)

