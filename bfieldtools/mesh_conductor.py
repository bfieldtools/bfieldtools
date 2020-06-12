"""
Contains a class used for wrapping a mesh (in the form of a Trimesh object) together with
some convenient functions and properties.

"""

__all__ = [
    "CouplingMatrix",
    "MeshConductor",
    "StreamFunction",
    "matrixwrapper",
    "save_pickle",
    "load_pickle",
]

from time import time
import pickle

import trimesh
import numpy as np
from scipy.sparse import issparse

from . import utils
from .mesh_calculus import laplacian_matrix, mass_matrix
from .mesh_impedance import (
    self_inductance_matrix,
    resistance_matrix,
    mutual_inductance_matrix,
)
from .mesh_magnetics import (
    magnetic_field_coupling,
    scalar_potential_coupling,
    vector_potential_coupling,
)
from .suhtools import SuhBasis
from .sphtools import compute_sphcoeffs_mesh
from .viz import plot_mesh, plot_data_on_vertices
from .contour import scalar_contour
from .line_conductor import LineConductor


def matrixwrapper(func):
    """
    Wrapper for lazy computation of MeshConductor matrices with basis change
    """

    def wrapper(obj):
        name = func.__name__
        M = obj.matrices[name[1:]]
        if M is None:
            print("Computing the %s matrix..." % str(name[1:]))
            M = obj.matrices[name[1:]] = func(obj)
        return obj.basis.T @ M @ obj.basis

    return wrapper


class MeshConductor:
    """
    Class that is used for surface mesh field calculations.
    
    Computation functions are typically external functions that are called
    using lazy properties.

    The mesh surface can consist of a single contiguous surface or several separate
    surfaces within a single mesh object. The MeshConductor object can handle data defined
    on the mesh being represented in several different bases:
     - inner (default)
     - vertex:
     - suh: surface harmonics basis. Order given by N_suh

    The bases can include built-in boundary conditions for the data: inner and
    suh bases assume dirichlet boundary condition (equal value within each boundary),
    while vertex basis does not set a boundary condition.
    
    """

    def __init__(
        self,
        verts=None,
        tris=None,
        mesh_file=None,
        mesh_obj=None,
        process=False,
        fix_normals=True,
        basis_name="inner",
        resistivity=1.68 * 1e-8,
        thickness=1e-4,
        **kwargs
    ):
        """
        First priority is to use given Trimesh object (mesh_obj).
        Second priority is to load mesh from file (mesh_file).
        Third priority is to use given verts and tris arrays (verts, tris).

        Parameters
        ----------
        verts: array-like (Nv, 3)
            Array of mesh vertices
        tris: array-like  (Nt, 3)
            Array of mesh faces
        mesh_obj: Trimesh mesh object
        mesh_file: string
            String describing the file path of a mesh file.
        process: boolean
            If True,  Trimesh will pre-process the mesh.
        fix_normals: boolean
            If True,  normals+winding should be set so that they always point "out" from the origin.
        basis_name: string
            Which basis to use, must be 'inner', 'vertex' or 'suh'. See class docstring
        Resistivity: float or array (Nfaces)
            Resistivity value in Ohm/meter
        Thickness: float or array (Nfaces)
            Thickness of surface. NB! Must be small in comparison to observation distance
        kwargs:
            Additional options with default settings are:
                'outer_boundaries':None,
                'resistance_full_rank': False,
                'inductance_nchunks':None,
                'basis_name':'inner' (other: suh, vertex)
                'N_suh': 100
                'sph_normalization': 'default'
                'sph_radius': 1
                'N_sph': 5
                'approx_far': True
                'approx_far_margin': 2
             
            
        Notes
        ------
        
        **outer_boundaries**
        int or array_like, indices of outer boundaries given by utils.find_boundaries(). 
        One boundary index per mesh component. If None, outer_boundaries are set to 
        the longest boundary in each mesh component. When using basis 'inner', the outer boundary
        vertex values are fixed to zero.
        
        **resistance_full_rank** (Boolean)
        If True, applies inflation to the resistance matrix in order to increase
        the rank by one. By default, is False. Don't set to True without knowing what you are doing.
        
        **inductance_nchunks** (int or None)
        Number of serial chunks to split self-inductance computation into, saving memory but taking more time.
        When approx_far is True, using more chunks is more efficient (multiply by 10-100x)
        If None (default), attempts to set number of chunks automatically based on the amount of free memory.
        Unfortunately, this estimation is not perfect.

        **N_suh**
        Number of surface harmonics to use if basis_name is 'suh'
        
        **sph_normalization**
        'default' (Ylm**2 integrated over solid angle to 1) or
        'energy' (field energy of basis fields normalized to 1 in R-ball)
        
        **sph_radius**
        If sph_normalization is 'energy', defines the radius of the inner expansion
        
        **N_sph**
        Number of spherical harmonics degrees (l-degrees) to use for the spherical harmonics coupling computation
        
        **approx_far** (Boolean)
        If True, usesimple quadrature for points far from the source triangles when computing self-inductance
        
        **approx_far_margin** (non-negative float)
        Cut-off distance for "far" points measured in mean triangle side length.


        """
        #######################################################################
        # Load mesh, do mesh pre-processing
        #######################################################################

        if (
            mesh_obj
        ):  # First, if Trimesh object is given as parameter then use that directly
            self.mesh = mesh_obj

        elif (
            mesh_file
        ):  # Second, check if mesh_file passed and load mesh as Trimesh object
            self.mesh = trimesh.load(mesh_file, process=process)

        else:  # Fall back on verts, tris parameters
            if isinstance(verts, type(None)) or isinstance(tris, type(None)):
                ValueError(
                    "You must provide either verts and tris, a mesh object or a mesh file"
                )

            self.mesh = trimesh.Trimesh(verts, tris, process=process)

        if fix_normals:
            self.mesh = utils.fix_normals(self.mesh)

        #######################################################################
        # Apply options
        #######################################################################

        # Populate options dictionary with defaults if not specified
        self.opts = {
            "outer_boundaries": None,
            "resistance_full_rank": False,
            "inductance_nchunks": None,
            "N_suh": 100,
            "sph_normalization": "default",
            "sph_radius": 1,
            "N_sph": 5,
            "inductance_quad_degree": 2,
            "approx_far": True,
            "approx_far_margin": 2,
        }

        for key, val in kwargs.items():
            self.opts[key] = val

        #######################################################################
        # Set up holes/boundaries
        #######################################################################

        self.boundaries, self.inner_vertices = utils.find_mesh_boundaries(self.mesh)
        self.set_holes(self.opts["outer_boundaries"])

        #######################################################################
        # Set up physical properties and coupling matrices
        #######################################################################

        # Matrices in inner-weight basis
        self.matrices = {
            "laplacian": None,
            "mass": None,
            "inductance": None,
            "resistance": None,
        }

        self.resistivity = resistivity
        self.thickness = thickness

        self.B_coupling = CouplingMatrix(self, magnetic_field_coupling)
        self.U_coupling = CouplingMatrix(self, scalar_potential_coupling)
        self.A_coupling = CouplingMatrix(self, vector_potential_coupling)

        # Coupling to spherical harmonic field coefficients
        self._alpha_coupling = None
        self._beta_coupling = None

        #######################################################################
        # Set up stream function basis
        #######################################################################

        self.inner2vert = utils.inner2vert(self.mesh, self.inner_vertices, self.holes)
        self.vert2inner = utils.vert2inner(self.mesh, self.inner_vertices, self.holes)

        # Sets basis for first time, calling self.set_basis()
        self.basis_name = basis_name
        self.set_basis(self.basis_name)

    def set_basis(self, basis_name):
        """ The data is stored in vertex basis i.e. every
            element corresponds to one vertex. The basis matrix changes basis
            of the operators so that a coefficient vector in the desired
            basis can be multiplied directly with the operator

            basis_names : str 'vertex', 'inner' or 'suh'
        """
        from scipy.sparse import spdiags

        self.basis_name = basis_name

        if self.basis_name == "suh":
            # SUH basis with dirichlet boundary conditions
            self.suh_basis = SuhBasis(
                self.mesh, self.opts["N_suh"], boundary_condition="dirichlet"
            )
            self.basis = self.inner2vert @ self.suh_basis.basis
        elif self.basis_name == "inner":
            self.basis = self.inner2vert
        elif self.basis_name == "vertex":
            N = len(self.mesh.vertices)
            self.basis = spdiags(np.ones(N), 0, N, N)
        else:
            raise ValueError("streamfunction_basis must inner, vertex or suh")

    def set_holes(self, outer_boundaries=None):
        """ Set indices of holes to self.holes

            outer_boundaries: int or array_like, indices of outer boundaries in
                    self.boundaries. One boundary index per mesh component.
                    If None, outer_boundaries are set the longest
                    boundary in each mesh component
        """
        if len(self.boundaries) == 0:
            # The mesh is watertight
            self.holes = []
            return
        if outer_boundaries is None:
            # Have to determine if there multiple bodies and label
            # the boundaries according to them
            comps = trimesh.graph.connected_components(self.mesh.edges)
            b_labels = np.zeros(len(self.boundaries))
            for m, b in enumerate(self.boundaries):
                for n, c in enumerate(comps):
                    if b[0] in c:
                        b_labels[m] = n
            b_lengths = np.array([len(b) for b in self.boundaries])
            # Determine outer_boundaries by finding the boundary of max length
            # for each component
            outer_boundaries = []
            b_inds = np.arange(len(self.boundaries))
            for n in range(len(comps)):
                mask = b_labels == n
                outer_boundaries.append(
                    b_inds[mask][np.argwhere(b_lengths[mask] == b_lengths[mask].max())]
                )

        self.opts["outer_boundaries"] = outer_boundaries
        hole_inds = list(
            np.setdiff1d(np.arange(len(self.boundaries)), outer_boundaries)
        )
        if len(hole_inds) == 0:
            # The non-watertight meshes contain no holes
            self.holes = []
        else:
            self.holes = [self.boundaries[i] for i in hole_inds]

    @property
    def laplacian(self):
        """
        Surface laplacian matrix, returned in appropiate basis.

        For further information, see mesh_calculus.laplacian_matrix

        property-decorated wrapper.
        """
        return self._laplacian()

    @matrixwrapper
    def _laplacian(self):
        """
        Compute and return surface laplacian matrix.

        """
        laplacian = laplacian_matrix(self.mesh, None)
        return laplacian

    @property
    def mass(self):
        """
        Mass matrix, returned in appropiate basis.

        For further information, see mesh_calculus.mass_matrix

        property-decorated wrapper.

        """
        return self._mass()

    @matrixwrapper
    def _mass(self):
        """
        Compute and return mesh mass matrix.

        """
        mass = mass_matrix(self.mesh, lumped=False)
        return mass

    @property
    def inductance(self):
        """
        Self-inductance matrix, returned in appropiate basis.

        For further information, see mesh_impedance.self_inductance_matrix

        property-decorated wrapper.
        """
        return self._inductance()

    @matrixwrapper
    def _inductance(self):
        """
        Compute and return self-inductance matrix.

        """

        start = time()

        inductance = self_inductance_matrix(
            self.mesh,
            Nchunks=self.opts["inductance_nchunks"],
            quad_degree=self.opts["inductance_quad_degree"],
            approx_far=self.opts["approx_far"],
            margin=self.opts["approx_far_margin"],
        )

        duration = time() - start
        print("Inductance matrix computation took %.2f seconds." % duration)

        return inductance

    @property
    def resistance(self):
        """
        Resistance matrix. For further information, see mesh_impedance.resistance_matrix

        property-decorated wrapper.
        """
        return self._resistance()

    @matrixwrapper
    def _resistance(self):
        """
        Back-end of resistance matrix computation
        """
        sheet_resistance = self.resistivity / self.thickness
        resistance = resistance_matrix(self.mesh, sheet_resistance).toarray()

        # Compensate for rank n-1 by adding offset, otherwise this
        # operator map constant vectors to zero
        if self.opts["resistance_full_rank"]:
            scale = np.mean(sheet_resistance)
            resistance += np.ones(resistance.shape) / resistance.shape[0] * scale

        return resistance

    def mutual_inductance(self, mesh_conductor_other, quad_degree=1, approx_far=True):
        """
        Mutual inductance between this MeshConductor object and another

        Parameters:
            mesh_conductor_other: MeshConductor object

        Returns:
            M: mutual inductance matrix M(self, other) in
                in the bases specified in the mesh_conductor object

        """
        M = mutual_inductance_matrix(
            self.mesh,
            mesh_conductor_other.mesh,
            quad_degree=quad_degree,
            approx_far=approx_far,
        )
        # Convert to the desired basis
        M = self.basis.T @ M @ mesh_conductor_other.basis

        return M

    def set_sph_options(self, **kwargs):
        """
        Set or reset options related to the spherical harmonics (SPH) coupling of the 
        MeshConductor object. If SPH couplings have already been computed, these will be 
        flushed.
    

        Parameters
        ----------
        **kwargs 
            Any parameter related to sph, namely "N_sph", "sph_normalization", "sph_radius" 

        Returns
        -------
        None.

        """

        # Flush old results
        self._alpha_coupling = None
        self._beta_coupling = None

        # Assign new values to opts, filter out opts not related to sph
        for key, val in kwargs.items():
            if key in ("N_sph", "sph_normalization", "sph_radius"):
                self.opts[key] = val

    @property
    def sph_couplings(self):
        """
        Spherical harmonic mappings from a stream function defined on the
        MeshConductor mesh.

        property-decorated wrapper
        """
        return self._sph_couplings()

    def _sph_couplings(self):
        """
        Compute spherical harmonic mappings and store them for further use
        """
        if self._alpha_coupling is None:
            print("Computing coupling matrices")
            Calpha, Cbeta = compute_sphcoeffs_mesh(
                self.mesh,
                self.opts["N_sph"],
                self.opts["sph_normalization"],
                self.opts["sph_radius"],
            )
            # Store the results for further use
            self._alpha_coupling = Calpha
            self._beta_coupling = Cbeta
        else:
            Calpha = self._alpha_coupling
            Cbeta = self._beta_coupling

        return Calpha @ self.basis, Cbeta @ self.basis

    def __setattr__(self, name, value):
        """
        Modified set-function to take into account post-hoc changes to e.g. resistance
        """
        self.__dict__[name] = value

        # If resistance-affecting parameter is changed after the resistance matrix has been computed,
        # then flush old result and re-compute
        if (name in ("resistivity", "thickness")) and self.matrices[
            "resistance"
        ] is not None:
            self.matrices["resistance"] = None  # Flush old matrix
            self._resistance()  # Re-compute with new parameters

        if (name in ("resistivity", "thickness")) and self.matrices[
            "resistance"
        ] is not None:
            self.matrices["resistance"] = None  # Flush old matrix
            self._resistance()  # Re-compute with new parameters

    def plot_mesh(self, cull_front=False, cull_back=False, **kwargs):
        """
        Simply plot the mesh surface in mayavi. kwargs are passed to
        viz.plot_mesh

        """

        return plot_mesh(
            self.mesh, cull_front=cull_front, cull_back=cull_back, **kwargs
        )

    def save_pickle(self, target_file):
        """
        Save the MeshConductor object using a pickled Python file

        Parameters
        ----------
        target_file: str
            File name or file object to save to

        """

        pickle.dump(obj=self, file=open(target_file, "wb"))


def save_pickle(obj, target_file):
    """
    Save the MeshConductor object using a pickled Python file

    Parameters
    ----------
    obj: object to save to file
    target_file: str
        file name or file object to save to

    """

    pickle.dump(obj=obj, file=open(target_file, "wb"), protocol=-1)


def load_pickle(target_file):
    """
    Load pickled MeshConductor object from file

    Parameters
    -----------
    target_file: str
        File name or file object to load from

    Returns
    -------
    obj: loaded MeshConductor object


    """

    obj = pickle.load(open(target_file, "rb"))

    return obj


class CouplingMatrix:
    """
    General-use class that contains a data array (a coupling matrix)
    and a bookkeeping list of computed points.

    When called, returns the coupling matrix for queried points.
    If some output has already been computed, use pre-computed values instead
    and only compute missing parts.


    """

    def __init__(self, parent, function):

        # Bookkeeping array, which points are already computed
        # Indexed in same order as the matrix
        self.points = np.array([])

        self.matrix = np.array([])

        self.parent = parent
        self.function = function

    def reset(self):
        """ Reset the coupling matrix and points
        """
        self.points = np.array([])
        self.matrix = np.array([])

    def __call__(self, points, s=None, *fun_args, **kwargs):
        """
        Returns the output of self.function(self.parent, points).
        If some output has already been computed, use pre-computed values instead
        and only compute missing parts.

        Parameters
        ----------
            points: (N_points, ... ) numpy array
                Array containing query points
            s: (Nbasis,) numpy array
                If None, return coupling matrix, 
                if not None, return field, i.e, coupling_matrix @ s
                default None
            *fun_args: additional arguments
                Optional, additional arguments that are passed to self.function

        Returns
        -------
            (N_points, ...) numpy array

        """

        if len(self.points) == 0:
            # Nothing calculated,yet
            self.matrix = self.function(self.parent.mesh, points, *fun_args, **kwargs)
            self.points = points

            M = self.matrix

        else:
            # Check which points exist and which need to be computed
            p_existing_point_idx, m_existing_point_idx = np.where(
                (self.points == points[:, None]).all(axis=-1)
            )
            missing_point_idx = np.setdiff1d(
                np.arange(0, len(points)), p_existing_point_idx
            )

            # If there are missing points, compute and add
            if len(missing_point_idx) > 0:
                missing_points = points[missing_point_idx]

                new_matrix_elems = self.function(
                    self.parent.mesh, missing_points, *fun_args, **kwargs
                )

                # Append newly computed point to coupling matrix, update bookkeeping
                self.points = np.vstack((self.points, missing_points))
                self.matrix = np.vstack((self.matrix, new_matrix_elems))

                # Re-compute indices of queried points, now that all should exist
                p_existing_point_idx, m_existing_point_idx = np.where(
                    (self.points == points[:, None]).all(axis=-1)
                )

            M = self.matrix[m_existing_point_idx]

        if s is None:
            if self.matrix.ndim == 2:
                M = M @ self.parent.basis

            elif self.matrix.ndim == 3:

                # Handle both sparse and dense basis matrices quickly
                if issparse(self.parent.basis):
                    Mnew = []
                    for n in range(3):
                        Mnew.append(M[:, n, :] @ self.parent.basis)
                    M = np.swapaxes(np.array(Mnew), 0, 1)
                else:
                    M = np.einsum("ijk,kl->ijl", M, self.parent.basis)

            else:
                raise ValueError("Matrix dimensions not ok")

            return M
        else:
            if self.matrix.ndim == 2:
                field = M @ (self.parent.basis @ s)
            elif self.matrix.ndim == 3:
                ss = self.parent.basis @ s
                field = np.einsum("ijk,k->ij", M, ss)
            else:
                raise ValueError("Matrix dimensions not ok")

            return field


class StreamFunction(np.ndarray):
    """ Class for representing stream function(s) on a MeshConductor

        Handles the mapping between different bases, e.g. inner vertices <->
        all vertices <-> surface harmonics

        Parameters:
            vals: array of shape (N,) or (N,M)
                where N corresponds to
                the number of inner vertices in the mesh_conductor or the
                the number of all vertices in the MeshConductor.

                Multiple (M) stream functions can be stored in the object
                by specifying vals with shape (N,M)
            mesh_conductor:
                MeshConductor object
    """

    def __new__(cls, input_array, mesh_conductor=None):
        #        print('In __new__ with class %s' % cls)
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.mesh_conductor = mesh_conductor

        return obj

    def __init__(self, input_array, mesh_conductor=None):
        #        print('In __init__ with class %s' % self.__class__)

        self.basis_name = self.mesh_conductor.basis_name
        self.basis = self.mesh_conductor.basis

        self.inner2vert = self.mesh_conductor.inner2vert
        self.vert2inner = self.mesh_conductor.vert2inner

    def __array_finalize__(self, obj):
        #
        #        print('In array_finalize:')
        #        print('   self type is %s' % type(self))
        #        print('   obj type is %s' % type(obj))

        if obj is None:
            #            obj = obj.view(np.ndarray)
            return

        self.mesh_conductor = getattr(obj, "mesh_conductor", None)

    def __array_wrap__(self, out_arr, context=None):
        """
        Return a ndarray for all ufuncs, since the StreamFunction attributes no longer apply
        when e.g. changing shape etc
        """
        return np.ndarray.__array_wrap__(self, out_arr, context).view(np.ndarray)

    def __getitem__(self, k):
        return np.ndarray.__getitem__(self, k).view(np.ndarray)

    def __repr__(self):
        return "%s, %s,  basis: %s" % (
            self.__class__.__name__,
            self.view(np.ndarray).__repr__(),
            self.basis_name,
        )

    @property
    def vert(self):
        """
        Returns the stream function in vertex basis
        """
        return self.basis @ self

    @property
    def inner(self):
        """
        Returns the stream function in inner basis
        """
        if self.basis_name == "inner":
            return self

        return self.vert2inner @ self.basis @ self

    @property
    def power(self):
        """
        Stream-function resistive power
        """

        # Compute resistance matrix if not present
        if self.mesh_conductor.matrices["resistance"] is None:
            self.mesh_conductor.resistance

        R = self.mesh_conductor.matrices["resistance"]  # Always in vertex basis

        return self.T @ self.basis.T @ R @ self.basis @ self

    @property
    def magnetic_energy(self):
        """
        Stream-function magnetic energy
        """

        # Compute inductance matrix if not present
        if self.mesh_conductor.matrices["inductance"] is None:
            self.mesh_conductor.inductance

        M = self.mesh_conductor.matrices["inductance"]  # Always in vertex basis

        return 0.5 * self.T @ self.basis.T @ M @ self.basis @ self

    def coil_inductance(self, Nloops):
        """
        

        Parameters
        ----------
        Nloops : the number of wire loops in the discrete coil (int)

        Returns
        -------
        inductance: scalar (float)

        """
        scaling = Nloops / (self.vert.max() - self.vert.min())
        L_approx = 2 * self.magnetic_energy * (scaling ** 2)

        return L_approx

    def plot(self, background=True, contours=False, **kwargs):
        """
        Plot the stream function
        """
        from mayavi import mlab

        mesh = self.mesh_conductor.mesh
        scalars = self.vert
        if "vmin" not in kwargs.keys():
            kwargs["vmin"] = -np.max(abs(scalars))
        if "vmax" not in kwargs.keys():
            kwargs["vmax"] = np.max(abs(scalars))

        s = plot_data_on_vertices(mesh, scalars, **kwargs)
        if contours:
            s.enable_contours = True
            s.contour.number_of_contours = contours
            if background:
                mlab.triangular_mesh(
                    *mesh.vertices.T, mesh.faces, color=(0.5, 0.5, 0.5), opacity=0.2
                )

        return s

    def discretize(self, N_contours=10, contours=None):
        """
        Wrapper method for contour.scalar_contour, turns the piecewise linear
        stream function into isolines/contours in the form of polylines.

        Parameters
        ----------
        N_contours: int
            Number of contours to generate
        contours: array-like
            Optional argument for manual input of contour levels. Overrides `N_contours`

        Returns
        -------
        contour_polys: list
            list with length `N_contours`. Each list element is anumpy array containing the
            coordinats of each polygon vertex.
        contour_values: array-like
            Vector containing the scalar function value for each contour line
        """
        return LineConductor(
            mesh=self.mesh_conductor.mesh,
            scalars=self.vert,
            N_contours=N_contours,
            contours=contours,
        )
