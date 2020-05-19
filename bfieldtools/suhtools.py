"""

Tools for calculating the surface harmonic representation of
the magnetic field. Surface harmonics can represent
of any stream function as a series.

Surface harmonics == Laplace-Beltrami eigenfunctions

"""

__all__ = ["SuhBasis"]

import trimesh
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh

import numpy as np

from .viz import plot_data_on_vertices
from . import mesh_conductor


class SuhBasis:
    """

    Class for representing magnetic field using surface harmonics

    """

    def __init__(
        self,
        obj,
        Nc=None,
        boundary_condition="dirichlet",
        magnetic=False,
        solver_sparse=True,
        **kwargs
    ):
        """
        Parameters
        ----------
        obj : Trimesh-object or Conductor-object
            Represents the boundary on which current density is specified
            or Conductor object that wraps the mesh
        Nc : Number of components
            If None (default), compute all components.
        boundary_condition : str  "dirichlet" (default) or "neumann"
            if zero-Dirichlet boundary conditions ("dirichlet")
            are used the basis corresponds to inner_vertices
            else with zero-Neumann condition ("neumann")
            the basis corresponds to all vertices
        magnetic: False or 'DC' or 'AC'
            Determines eigenvalue equation. If False, use laplacian and mass matrices.
            If 'DC', use resistance and mass matrices. If 'AC', use resistance and inductance matrices.
        solver_sparse: Boolean (True)
            If True, use solver from scipy.sparse.linalg rather than scipy.linalg
        kwargs: dict
            Passed to Conductor creation if a Trimesh object is passed as 'obj'

        """

        if boundary_condition in ("dirichlet", "neumann"):
            self.bc = boundary_condition
        else:
            raise ValueError(
                "boundary_conditions should be either dirichlet or neumann"
            )

        if isinstance(obj, mesh_conductor.MeshConductor):
            if obj.opts["resistance_full_rank"] is True:
                raise ValueError(
                    "resistance matrix must not be deflated (resistance_full_rank must be False)"
                )
            self.mesh_conductor = obj
        elif isinstance(obj, trimesh.Trimesh):
            self.mesh_conductor = mesh_conductor.MeshConductor(
                mesh_obj=obj, resistance_full_rank=False, **kwargs
            )
        else:
            raise TypeError("obj type should be either Trimesh or Conductor")
        self.mesh = self.mesh_conductor.mesh

        self.magnetic = magnetic
        self.solver_sparse = solver_sparse

        if self.bc == "neumann":
            self.mesh_conductor.set_basis("vertex")
        else:
            self.mesh_conductor.set_basis("inner")

        if self.mesh.is_watertight or self.bc == "neumann":
            self._max_Nc = self.mesh_conductor.basis.shape[1] - 1
        else:
            self._max_Nc = self.mesh_conductor.basis.shape[1]

        if Nc is None:
            print("No component count given, computing all components.")
            self.Nc = self._max_Nc
        elif Nc > self._max_Nc:
            raise ValueError("Mesh supports up to %d components" % self._max_Nc)
        else:
            self.Nc = Nc

        self.inner_vertices = self.mesh_conductor.inner_vertices
        self.holes = self.mesh_conductor.holes
        self.inner2vert = self.mesh_conductor.inner2vert

        self.calculate_basis()

    def calculate_basis(self, shiftinvert=True, v0=None):
        """ Calculate basis functions as eigenfunctions of matrices
        specified by the 'magnetic' parameter
        
        Parameters
        ----------
        shiftinvert: Boolean (True)
            use shiftinvert mode to calculate eigenstuff faster (experimental)
        v0: array (b, b)
            Intial value matrix for eigenvalue decomposition algorithm. Shape depends on BC: if neumann b=self.mesh_conductor.basis.shape[0],
            if dirichlet b=self.mesh_conductor.basis.shape[1]
                        
                    
        """
        print("Calculating surface harmonics expansion...")

        if not self.magnetic:
            L = -self.mesh_conductor.laplacian
            M = self.mesh_conductor.mass
        elif self.magnetic == "AC":
            L = self.mesh_conductor.resistance
            M = self.mesh_conductor.inductance
            # If mesh is closed, add 'deflation' to the inductance
            # L matrix deflation is already provided by Conductor(resistance_full_rank=True)
            if self.mesh.is_watertight:
                M += np.ones(M.shape) * np.mean(np.diag(M))
        elif self.magnetic == "DC":
            L = self.mesh_conductor.resistance
            M = self.mesh_conductor.mass
        else:
            raise ValueError("Parameter 'magnetic' must be False, 'DC', or 'AC'.")

        if self.mesh.is_watertight or self.bc == "neumann":
            print("Closed mesh or Neumann BC, leaving out the constant component")

            N0 = 1

            N = np.min((self.Nc + 1, self._max_Nc))

            # adjust self.Nc so that is corresponds to suh.basis.shape[1]
            if N == self._max_Nc:
                self.Nc -= 1
        else:
            N0 = 0
            N = self.Nc

        # Make sure that matrices are symmetric, positive definite
        L = 0.5 * (L + L.T)
        M = 0.5 * (M + M.T)

        if v0 is None:
            v0 = np.ones(L.shape[1])  # avoid random basis for symmetric geometries

        # Use sparse solver if possible
        # 1. matrices have to be sparse (i.e. magnetic=False)
        # 2. solver_sparse = True
        # 3. the sparse solver cannot find all eigenvalues, so Nc < max_Nc
        if (not self.magnetic) and self.solver_sparse and (self.Nc < self._max_Nc):

            if shiftinvert:
                u, v = eigsh(L, N, M, sigma=0, which="LM", v0=v0)
            else:
                u, v = eigsh(L, N, M, which="SA", v0=v0)

            self.basis = v[:, N0:]
            self.eigenvals = u[N0:]

        else:
            if not self.magnetic:
                L = L.toarray()
                M = M.toarray()

                u, v = eigh(L, M, eigvals=(0, N - 1))
            else:
                if self.magnetic == "DC":
                    M = M.toarray()

                u, v = eigh(L, M, eigvals=(0, N - 1))

            self.basis = v[:, N0:]
            self.eigenvals = u[N0:]

    def field(self, coeffs, points):
        """ Calculate field at points

            Parameters

            coeffs : (self.Nc,) array of basis function coefficients
            points : (N_points, 3) field evaluation points

            Returns:

                field : (N_points, 3) magnetic field
        """
        return self.basis_fields(points) @ coeffs

    def basis_fields(self, points):
        """ Calculate basis fields at points

            Return:

                Fields (3, N_points, self.Nc)
        """
        B_coupling = self.mesh_conductor.B_coupling(points)

        return B_coupling @ self.basis

    def fit_coeffs(self, points, data):
        """ Fit basis function coefficients to the data
        """
        assert len(data) > self.Nc

        A = self.basis_fields(points).reshape(-1, self.Nc)
        b = data.T.flatten()
        x, res, rank, s = np.linalg.lstsq(A, b, rcond=None)

        if rank < self.Nc:
            print("Matrix rank not full, result might be inaccurate")
        return x

    def plot(
        self, Nfuncs, dist=0.5, Ncols=None, figure=None, figsize=(800, 800), **kwargs
    ):
        """ Plot basis functions on the mesh

            Nfuncs: int or array-like
                 if int, the number functions starting from the first'Blm',
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
        from mayavi import mlab

        if figure is None:
            figure = mlab.figure(
                None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=figsize
            )

        if isinstance(Nfuncs, int):
            N = Nfuncs
            indices = np.arange(Nfuncs)
        else:
            indices = Nfuncs
            N = len(indices)

        if Ncols is None:
            Ncols = np.floor(np.sqrt(N) + 1)

        dx = (self.mesh.vertices[:, 0].max() - self.mesh.vertices[:, 0].min()) * (
            1 + dist
        )
        dy = (self.mesh.vertices[:, 1].max() - self.mesh.vertices[:, 1].min()) * (
            1 + dist
        )

        i = 0
        j = 0

        s = []
        for n in indices:
            print(i, j)

            tmp_mesh = self.mesh.copy()
            tmp_mesh.vertices[:, 0] += i * dx
            tmp_mesh.vertices[:, 1] -= j * dy

            if self.bc == "neumann":
                scalars = self.basis[:, n]
            else:
                scalars = self.inner2vert @ self.basis[:, n]

            s.append(plot_data_on_vertices(tmp_mesh, scalars, figure=figure, **kwargs))

            if i < Ncols:
                i += 1
            else:
                j += 1
                i = 0

        return s
