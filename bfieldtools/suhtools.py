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

from .mesh_calculus import laplacian_matrix, mass_matrix
from .mesh_magnetics import magnetic_field_coupling
from .utils import inner2vert
from scipy.sparse.linalg import eigsh
import numpy as np
from mayavi import mlab

class SuhBasis():
    """ Class for representing magnetic field using surface harmonics
    """

    def __init__(self, mesh, Nc, closed_mesh=True, inner_vertices=None, holes=None):
        """
        Parameters
            mesh : Trimesh-object representing the boundary on which
                        current density is specified
            Nc : Number of components
            closed_mesh:
                inform if the mesh is closed or not (could be automated)
            inner_vertices
                If given zero-Dirichlet boundary conditions are used
                for the calculation, other all vertices are used
                which corresponds to the (natural) zero-Neumann condition
            holes:
                list of lists of indices for vertices that belong to holes
        """

        self.mesh = mesh
        self.Nc = Nc
        if inner_vertices is not None:
            self.inner_vertices = inner_vertices
        else:
            self.inner_vertices = np.arange(len(self.mesh.vertices))
        self.holes = holes
        self.calculate_basis(closed_mesh)
        if self.holes is not None:
            self.inner2vert = inner2vert(self.mesh, self.inner_vertices, self.holes)
        elif inner_vertices is not None:
            self.inner2vert = inner2vert(self.mesh, self.inner_vertices, [])
            self.holes = []
        else:
            self.inner2vert = np.eye(len(self.mesh.vertices))

    def calculate_basis(self, closed_mesh=True, shiftinvert=False):
        """ Calculate basis functions as eigenfunctions of the laplacian

            closed_mesh: if True, calculate the basis for the whole mesh
                         if False, calculate for inner vertices (zero-dirichlet condition)
                                    and use zero for the boundary
        """
        print('Calculating surface harmonics expansion...')

        if closed_mesh:
            assert self.mesh.is_watertight

        L = laplacian_matrix(self.mesh, None, self.inner_vertices, self.holes)
        M = mass_matrix(self.mesh, False, self.inner_vertices, self.holes)
        self.mass = M

        if closed_mesh:
            N0 = 1
            N = self.Nc + 1
        else:
            N0 = 0
            N = self.Nc

        v0 = np.ones(L.shape[1]) # avoid random basis for symmetric geometries
        if shiftinvert:
            u, v = eigsh(-L, N, M, sigma=0, which='LM', v0 = v0)
        else:
            u, v = eigsh(-L, N, M, which='SA', v0 = v0)

        # The first function is constant and does not yield any field
        self.basis = v[:,N0:]
        self.eigenvals = u[N0:]

    def field(self, coeffs, points):
        """ Calculate field at points

            Parameters

            coeffs : (self.Nc,) array of basis function coefficients
            points : (N_points, 3) field evaluation points

            Returns:

                field : (N_points, 3) magnetic field
        """
        return (self.basis_fields(points) @ coeffs)

    def basis_fields(self, points):
        """ Calculate basis fields at points

            Return:

                Fields (3, N_points, self.Nc)
        """
        B_coupling = magnetic_field_coupling(mesh, points)

        return B_coupling @ self.basis

    def fit_coeffs(self, points, data):
        """ Fit basis function coefficients to the data
        """
        assert len(data) > self.Nc

        A = self.basis_fields(points).reshape(-1, self.Nc)
        b = data.T.flatten()
        x, res, rank, s = np.linalg.lstsq(A, b, rcond=None)

        if rank < self.Nc:
            print('Matrix rank not full, result might be inaccurate')
        return  x

    def plot(self, Nfuncs, dist=0.5, ncol=None, ncolors=15, **kwargs):
        """ Plot basis functions on the mesh

            Nfuncs: int or array-like
                 if int, the number functions starting from the first
                 if list/array: the indices of the functions

            dist: float
                distance between the plotted objects relative to their size

            ncol: int or None
                the number of columns in the plot
                If none automatically determined

            ncolors:
                number of colors in the colormap

            kwargs: keyword arguments passed to mayavi (colormap, etc.)

        """
        if type(Nfuncs) == int:
            N=Nfuncs
            indices = np.arange(Nfuncs)
        else:
            indices = Nfuncs
            N = len(indices)

        if ncol is None:
            ncol = np.floor(np.sqrt(N)+1)
        dx = (self.mesh.vertices[:,0].max() - self.mesh.vertices[:,0].min())*(1+dist)
        dy = (self.mesh.vertices[:,1].max() - self.mesh.vertices[:,1].min())*(1+dist)

        i = 0
        j = 0

        for n in indices:
            print(i,j)
            points = self.mesh.vertices.copy()
            points[:,0] += i*dx
            points[:,1] -= j*dy
            scalars = self.inner2vert @ self.basis[:,n]
            s = mlab.triangular_mesh(*points.T, self.mesh.faces,
                                 scalars=scalars, **kwargs)

            s.module_manager.scalar_lut_manager.number_of_colors = ncolors
            s.actor.mapper.interpolate_scalars_before_mapping = True
            if i<ncol:
                i+=1
            else:
                j+=1
                i=0

        return s


class SuhBasis2(SuhBasis):
    '''
    Implement more general version here
    '''
    def __init__(self):
        pass

if __name__ == '__main__':
    """ Simple testing script
    """

    from trimesh.creation import icosphere

    # Create basis for a sphere (basis.eigenvals shere the same structure
    # as spherical harmonic eigenvalues)
#    mesh = icosphere(4)
    import pkg_resources
    import trimesh
    #Load simple plane mesh that is centered on the origin
    file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/closed_cylinder_remeshed.stl')
    mesh = trimesh.load(file_obj, process=True)
    basis = SuhBasis(mesh, 40, True)

#    s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces,
#                             scalars=basis.basis[:,6])
#    s.enable_contours = True
#
#    coeffs = np.zeros(basis.basis.shape[1])
#    coeffs[6] = 1
#
#    # Plot field outside
#    points = mesh.vertices * 2.0
#    field = basis.field(coeffs, points)
#    v = mlab.quiver3d(*points.T, *field.T)
#    v.glyph.glyph_source.glyph_position = 'center'
#
#    # Plot field inside
#    points = mesh.vertices * 0.5
#    field = basis.field(coeffs, points)
#    v = mlab.quiver3d(*points.T, *field.T)
#    v.glyph.glyph_source.glyph_position = 'center'
#
#    # Modify the field and fit coeffs
#    field1 = np.zeros_like(field)
#    field1[:, 1] = np.max(field)
#    field2 = field + field1
#
#    coeffs2 = basis.fit_coeffs(points, field2)
#
#    mlab.figure()
#    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=basis.basis @ coeffs2)

