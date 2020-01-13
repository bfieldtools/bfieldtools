'''
Contains a class used for wrapping a mesh (in the form of a Trimesh object) together with
some convenient functions and properties.

'''

from time import time
import pickle

from mayavi import mlab
import trimesh
import numpy as np

from . import utils
from .mesh_calculus import laplacian_matrix, mass_matrix
from .mesh_properties import self_inductance_matrix, resistance_matrix
from .mesh_magnetics import magnetic_field_coupling, scalar_potential_coupling, vector_potential_coupling
from .suhtools import SuhBasis


class LazyProperty():
    '''
    Implementation of lazily loading properties, see
    http://blog.pythonisito.com/2008/08/lazy-descriptors.html
    On first invocation, a lazy property calls a function that populates
    the property (acts as a method). Afterwards, it acts like a normal property.

    '''

    def __init__(self, func):
        self._func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__

    def __get__(self, obj, klass=None):
        if obj is None:
            return None

        result = obj.__dict__[self.__name__] = self._func(obj)

        return result


class Conductor:
    '''
    Class that is used for surface mesh field calculations, e.g. coil design.
    Computation functions are typically external functions that are called
    using lazy properties.

    The mesh surface can consist of a single contiguous surface or several separate
    surfaces within a single mesh object.

    '''

    def __init__(self, verts=None, tris=None, mesh_file=None,
                 mesh_obj=None, process=False, fix_normals=True,
                 resistivity=1.68*1e-8, thickness=1e-4,
                 opts = {}):
        '''
        Initialize Conductor object.
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
        Resistivity: float or array (Nfaces)
            Resistivity value in Ohm/meter
        Thickness: float or array (Nfaces)
            Thickness of surface. NB! Must be small in comparison to observation distance
        opts: dict
            Options for Conductor object. Default settings are:
                'outer_boundaries':None,
                'mass_lumped':False,
                'resistance_full_rank': True,
                'inductance_nchunks':None,
                'basis':'free' (other: suh, vertex)


        '''
        #######################################################################
        # Load mesh, do mesh pre-processing
        #######################################################################

        if mesh_obj: #First, if Trimesh object is given as parameter then use that directly
            self.mesh = mesh_obj

        elif mesh_file: #Second, check if mesh_file passed and load mesh as Trimesh object
            self.mesh = trimesh.load(mesh_file, process=process)

        else: #Fall back on verts, tris parameters
            if isinstance(verts, type(None)) or isinstance(tris, type(None)):
                ValueError('You must provide either verts and tris, a mesh object or a mesh file')

            self.mesh = trimesh.Trimesh(verts, tris, process=process)

        if fix_normals:
            self.mesh = utils.fix_normals(self.mesh)

        #######################################################################
        # Apply options
        #######################################################################

        #Populate options dictionary with defaults if not specified
        self.opts = {'outer_boundaries':None, 'mass_lumped':False,
                     'resistance_full_rank': True, 'inductance_nchunks':None,
                     'streamfunction_basis':'vertex'}

        for key, val in opts.items():
            self.opts[key] = val


        #######################################################################
        # Set up holes/boundaries
        #######################################################################

        self.boundaries, self.inner_verts = utils.find_mesh_boundaries(self.mesh)
        self.set_holes(self.opts['outer_boundaries'])

        #######################################################################
        # Set up stream function basis
        #######################################################################

        self.f2v = utils.free2vert(self.mesh, self.inner_verts, self.holes)
        self.v2f = self.f2v.T

        self.basis_name = self.opts['streamfunction_basis']
        self.set_basis()

        #######################################################################
        # Set up physical properties and coupling matrices
        #######################################################################


        self.__dict__['resistivity'] = resistivity
        self.__dict__['thickness'] = thickness

        self.B_coupling = CouplingMatrix(self, magnetic_field_coupling)
        self.U_coupling = CouplingMatrix(self, scalar_potential_coupling)
        self.A_coupling = CouplingMatrix(self, vector_potential_coupling)


    def set_basis(self):
        if self.opts['streamfunction_basis'] == 'suh':
            self.basis = SuhBasis(self.mesh, self.mesh.is_watertight,
                                  self.inner_vertices, self.holes)
        elif self.opts['streamfunction_basis'] == 'free':
            self.basis = np.eye(len(self.inner_verts) + len(self.holes()))
        elif self.opts['streamfunction_basis'] == 'vertex':
            self.basis = self.v2f
        else:
            raise ValueError('streamfunction_basis must free, vertex or suh')

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
            # Determine outer_boundaries by finding the boundary of max lenght
            # for each component
            outer_boundaries = []
            b_inds = np.arange(len(self.boundaries))
            for n in range(len(comps)):
                mask = b_labels == n
                outer_boundaries.append(b_inds[mask][b_lengths[mask].argmax()])

        self.opts['outer_boundaries'] = outer_boundaries
        hole_inds = list(np.setdiff1d(np.arange(len(self.boundaries)),
                                 outer_boundaries))
        if len(hole_inds) == 0:
            # The non-watertight meshes contain no holes
            self.holes = []
        else:
            self.holes = self.boundaries[hole_inds]


    @LazyProperty
    def laplacian(self):
        '''
        Compute and return surface laplacian matrix.

        '''
        if len(self.holes) == 0:
            laplacian = laplacian_matrix(self.mesh, None, self.inner_verts)
        else:
            laplacian = laplacian_matrix(self.mesh, None, self.inner_verts,
                                         self.holes)
        return laplacian


    @LazyProperty
    def mass(self):
        '''
        Compute and return mesh mass matrix.

        '''
        if len(self.holes) == 0:
            mass = mass_matrix(self.mesh, self.opts['mass_lumped'], self.inner_verts)
        else:
            mass = mass_matrix(self.mesh, self.opts['mass_lumped'],
                               self.inner_verts, self.holes)

        return mass


    @LazyProperty
    def inductance(self):
        '''
        Compute and return mutual inductance matrix.

        '''

        start = time()

        inductance = self_inductance_matrix(self.mesh,
                                            Nchunks=self.opts['inductance_nchunks'])

        duration = time() - start
        print('Inductance matrix computation took %.2f seconds.'%duration)

        U = self.v2f
        return U @ inductance @ U.T


    def __setattr__(self, name, value):
        '''
        Modified set-function to take into account post-hoc changes to e.g. resistance
        '''
        self.__dict__[name] = value

        #If resistance-affecting parameter is changed after the resistance matrix has been computed,
        #then flush old result and re-compute
        if (name == "resistivity" or name == "thickness") and 'resistance' in self.__dict__.keys():
            self.resistance = self._resistance() #Re-compute with new parameters

        if name == 'basis_name':
            self.set_basis()


    def __getattr_(self, name):
        '''
        Modified get-function to implement basis mapping
        '''

        if name in ('laplacian', 'mass', 'resistance', 'inductance'):
            M = self.__dict__[name]
            U = self.basis
            return U.T @ M @ U

        return self.__dict__[name]


    @LazyProperty
    def resistance(self):
        '''
        Compute and return resistance/resistivity matrix using Laplace matrix.
        Conductivity and thickness are class parameters

        Returns
        -------
        R: array
            Resistance matrix

        '''

        return self._resistance()


    def _resistance(self):
        '''
        Back-end of resistance matrix computation
        '''
        sheet_resistance = self.resistivity / self.thickness
        resistance =  resistance_matrix(self.mesh, sheet_resistance).toarray()

        # Compensate for rank n-1 by adding offset, otherwise this
        # operator map constant vectors to zero
        if self.opts['resistance_full_rank']:
            scale = np.mean(sheet_resistance)
            resistance += np.ones(resistance.shape)/resistance.shape[0]*scale

        U = self.v2f
        return U @ resistance @ U.T



    def plot_mesh(self, representation='wireframe', opacity=0.5, color=(0, 0, 0), cull_front=False, cull_back=False):
        '''
        Simply plot the mesh surface in mayavi.

        '''

        mesh = mlab.triangular_mesh(*self.mesh.vertices.T, self.mesh.faces,
                                    representation=representation, opacity=opacity, color=color)

        mesh.actor.property.frontface_culling = cull_front
        mesh.actor.property.backface_culling = cull_back
        return mesh


    def save_pickle(self, target_file):
        """
        Save the Conductor object using a pickled Python file

        Parameters
        ----------
        target_file: str
            File name or file object to save to

        """

        pickle.dump(obj=self, file=open(target_file, 'wb'))


def save_pickle(obj, target_file):
    """
    Save the Conductor object using a pickled Python file

    Parameters
    ----------
    obj: object to save to file
    target_file: str
        file name or file object to save to

    """

    pickle.dump(obj=obj, file=open(target_file, 'wb'), protocol=-1)


def load_pickle(target_file):
    """
    Load pickled Conductor object from file

    Parameters
    -----------
    target_file: str
        File name or file object to load from
    Returns
    -------
    obj: loaded Conductor object


    """

    obj = pickle.load(open(target_file, 'rb'))

    return obj


class CouplingMatrix:
    '''
    General-use class that contains a data array (a coupling matrix)
    and a bookkeeping list of computed points.

    When called, returns the coupling matrix for queried points.
    If some output has already been computed, use pre-computed values instead
    and only compute missing parts.


    '''
    def __init__(self, parent, function):

        #Bookkeeping array, which points are already computed
        #Indexed in same order as the matrix
        self.points = np.array([])

        self.matrix = np.array([])

        self.parent = parent
        self.function = function


    def __call__(self, points, *fun_args):
        '''
        Returns the output of self.function(self.parent, points).
        If some output has already been computed, use pre-computed values instead
        and only compute missing parts.

        Parameters
        ----------
            points: (N_points, ... ) numpy array
                Array containing query points
            *fun_args: additional arguments
                Optional, additional arguments that are passed to self.function

        Returns
        -------
            (N_points, ...) numpy array

        '''

        if len(self.points) == 0:
            self.matrix = self.function(self.parent.mesh, points, *fun_args)
            # Convert to all-vertices to free vertices
            self.matrix = self.matrix @ self.parent.v2f
            self.points = points

            M = self.matrix

        else:
            #Check which points exist and which need to be computed
            p_existing_point_idx, m_existing_point_idx = np.where((self.points == points[:, None]).all(axis=-1))
            missing_point_idx = np.setdiff1d(np.arange(0, len(points)), p_existing_point_idx)

            #If there are missing points, compute and add
            if len(missing_point_idx) > 0:
                missing_points = points[missing_point_idx]

                new_matrix_elems = self.function(self.parent.mesh, missing_points, *fun_args)
                new_matrix_elems = new_matrix_elems @ self.parent.v2f


                #Append newly computed point to coupling matrix, update bookkeeping
                self.points = np.vstack((self.points, missing_points))
                self.matrix = np.vstack((self.matrix, new_matrix_elems))

                #Re-compute indices of queried points, now that all should exist
                p_existing_point_idx, m_existing_point_idx = np.where((self.points == points[:, None]).all(axis=-1))

            M = self.matrix[m_existing_point_idx]

        return M @ self.parent.basis


class StreamFunction:
    """ Class for representing stream function(s) on a conductor

        Handles the mapping between the degrees of freedom in the
        stream function (dof) and the vertex weights (w)

        Parameters:
            vals:
                    array of shape (N,) or (N,M) where N corresponds to
                    the number of free vertices in the conductor or the
                    the number of all vertices in the conductor.

                Multiple (M) stream functions can be stored in the object
                by specifying vals with shape (N,M)
            conductor:
                Conductor object
    """
    def __init__(self, vals, conductor):
        self.conductor = conductor
        self.f2v = self.conductor.f2v
        self.v2f = self.conductor.v2f
        self.set_stream_func(vals)

    def set_stream_func(self, vals):
        """ Set stream function values to the object

            Can also be used for re-setting the values

            Parameters:
                vals:
                    array of shape (N,) or (N,M) where N corresponds to
                    the number of free vertices in the conductor or the
                    the number of all vertices in the conductor.
        """
        if len(vals) == len(self.conductor.mesh.vertices):
            self.free = self.v2f @ vals
        elif len(vals) == len(self.conductor.inner_vertices) + len(self.conductor.holes):
            self.free = vals
        else:
            raise ValueError('the length of vals must either correspond to that of free vertices or vertex weights')


    @property
    def w(self):
        return self.v2f @ self.free

    @property
    def f(self):
        return self.free

    @property
    def power(self):
        return self.d.T @ self.conductor.resistance @ self.d

    @property
    def magnetic_energy(self):
        return 0.5 * self.d.T @ self.conductor.inductance @ self.d

    def plot(self):
        pass








