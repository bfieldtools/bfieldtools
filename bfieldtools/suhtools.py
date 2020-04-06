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
from .viz import plot_data_on_vertices
from . import conductor
import trimesh 

from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh

import numpy as np
from mayavi import mlab

class SuhBasis():
    """ Class for representing magnetic field using surface harmonics
    """

    def __init__(self, obj, Nc, boundary_condition='dirichlet',
                 magnetic=False, solver_sparse=True):
        """
        Parameters
            obj : Trimesh-object representing the boundary on which
                        current density is specified
                  or Conductor object that wraps the mesh
            Nc : Number of components
            boundary_condition : str  "dirichlet" (default) or "neumann"
                if zero-Dirichlet boundary conditions ("dirichlet") 
                are used the basis corresponds to inner_vertices
                else with zero-Neumann condition ("neumann") 
                the basis corresponds to all vertices
              
        """
        
        if boundary_condition in ("dirichlet", "neumann"):
          self.bc = boundary_condition
        else:
          raise ValueError('boundary_conditions should e either dirichlet or neumann')

        if isinstance(obj, conductor.Conductor):
          self.conductor = obj
        elif isinstance(obj, trimesh.Trimesh):
          # TODO defaults ok?
          self.conductor = conductor.Conductor(mesh_obj=obj, resistance_full_rank=False)
        else:
          raise TypeError('obj type should be either Trimesh or Conductor')
        self.mesh = self.conductor.mesh

        self.Nc = Nc
        self.magnetic = magnetic
        self.solver_sparse = solver_sparse
        
        
        if self.bc == 'neumann':
          self.conductor.set_basis('vertex')
        else:
          self.conductor.set_basis('inner')

          self.inner_vertices = self.conductor.inner_vertices
          self.holes = self.conductor.holes
          self.inner2vert = self.conductor.inner2vert
          
        
        self.calculate_basis()


    def calculate_basis(self, shiftinvert=True, v0=None):
        """ Calculate basis functions as eigenfunctions of the laplacian

            shiftinvert: use shiftinvert mode to calculate eigenstuff faster
                        (experimental)
        """
        print('Calculating surface harmonics expansion...')


        if not self.magnetic:
          L = self.conductor.laplacian
          M = self.conductor.mass
        else:
          L = -self.conductor.resistance
          M = self.conductor.inductance
          
        self.mass = M

        closed_mesh = self.conductor.mesh.is_watertight
        
        if closed_mesh:
            print('Closed mesh, leaving out the constant component')
            N0 = 1
            N = self.Nc + 1
        else:
            N0 = 0
            N = self.Nc
        
        if v0 is None:
            v0 = np.ones(L.shape[1]) # avoid random basis for symmetric geometries
            
        if (not self.magnetic) and self.solver_sparse:         
          if shiftinvert:
              u, v = eigsh(-L, N, M, sigma=0, which='LM', v0 = v0)
          else:
              u, v = eigsh(-L, N, M, which='SA', v0 = v0)
        else:
            if not self.magnetic:
                L = L.toarray()
                M = M.toarray()
            u, v = eigh(-L, M, eigvals=(0, N))

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


    def plot(self, Nfuncs, dist=0.5, Ncols=None, figure=None,
             figsize=(800,800), **kwargs):
        """ Plot basis functions on the mesh

            Nfuncs: int or array-like
                 if int, the number functions starting from the first
                 if list/array: the indices of the functions

            dist: float
                distance between the plotted objects relative to their size

            Ncols: int or None
                the number of columns in the plot
                If none automatically determined

            fig: handle for mlab figure

            figsize: size of a new figure if 'fig' not given

            ncolors:
                number of colors in the colormap

            kwargs: keyword arguments passed to mayavi (colormap, etc.)

        """

        if figure is None:
            figure = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                             size=figsize)

        if type(Nfuncs) == int:
            N=Nfuncs
            indices = np.arange(Nfuncs)
        else:
            indices = Nfuncs
            N = len(indices)

        if Ncols is None:
            Ncols = np.floor(np.sqrt(N)+1)

        dx = (self.mesh.vertices[:,0].max() - self.mesh.vertices[:,0].min())*(1+dist)
        dy = (self.mesh.vertices[:,1].max() - self.mesh.vertices[:,1].min())*(1+dist)

        i = 0
        j = 0

        for n in indices:
            print(i,j)

            tmp_mesh = self.mesh.copy()
            tmp_mesh.vertices[:,0] += i*dx
            tmp_mesh.vertices[:,1] -= j*dy

            scalars = self.inner2vert @ self.basis[:,n]

            s = plot_data_on_vertices(tmp_mesh, scalars, figure=figure, **kwargs)

            if i<Ncols:
                i+=1
            else:
                j+=1
                i=0

        return s



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
    basis = SuhBasis(mesh, 40)

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

