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

        '''

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

        self.opts = {'outer_boundaries':None, 'mass_lumped':False,
                     'resistance_full_rank': True, 'outer_boundaries':None}
        for key, val in opts.items():
            self.opts[key] = val


        self.boundaries, self.inner_verts = utils.find_mesh_boundaries(self.mesh)

        self.set_holes(self.opts['outer_boundaries'])

        self.resistivity = resistivity
        self.thickness = thickness

        self.B_coupling = CouplingMatrix(self, magnetic_field_coupling)
        self.U_coupling = CouplingMatrix(self, scalar_potential_coupling)
        self.A_coupling = CouplingMatrix(self, vector_potential_coupling)


        self.s = None
        self.problem = None


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
            laplacian = laplacian_matrix(self.mesh)
        else:
            laplacian = laplacian_matrix(self.mesh, None, self.inner_vertices,
                                         self.holes)
        return laplacian


    @LazyProperty
    def mass(self):
        '''
        Compute and return mesh mass matrix.

        '''
        if len(self.holes) == 0:
            mass = mass_matrix(self.mesh, self.opt['mass_lumped'])
        else:
            mass = mass_matrix(self.mesh, self.opt['mass_lumped'],
                               self.inner_vertices, self.holes)

        return mass


    @LazyProperty
    def inductance(self):
        '''
        Compute and return mutual inductance matrix.

        '''

        start = time()

        # TODO opts?
        inductance = self_inductance_matrix(self.mesh)# , Nchunks=self.opts['Nchunks'])

        duration = time() - start
        print('Inductance matrix computation took %.2f seconds.'%duration)

        return inductance

    @LazyProperty
    def resistance(self):
        '''
        Compute and return resistance/resistivity matrix using Laplace matrix.
        Default resistivity set to that of copper.
        Parameters
        ----------
        Resistivity: float or array (Nfaces)
            Resistivity value in Ohm/meter
        Thickness: float or array (Nfaces)
            Thickness of surface. NB! Must be small in comparison to observation distance

        Returns
        -------
        R: array
            Resistance matrix

        '''

        sheet_resistance = self.resistivity / self.thickness
        resistance =  resistance_matrix(self.mesh, sheet_resistance).todense()

        # Compensate for rank n-1 by adding offset, otherwise this
        # operator map constant vectors to zero
        if self.opts['resistance_full_rank']:
            scale = np.mean(sheet_resistance)
            resistance += np.ones(resistance.shape)/resistance.shape[0]*scale

        return resistance


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
            self.points = points

            return self.matrix

        else:
            #Check which points exist and which need to be computed
            p_existing_point_idx, m_existing_point_idx = np.where((self.points == points[:, None]).all(axis=-1))
            missing_point_idx = np.setdiff1d(np.arange(0, len(points)), p_existing_point_idx)

            #If there are missing points, compute and add
            if len(missing_point_idx) > 0:
                missing_points = points[missing_point_idx]

                new_matrix_elems = self.function(self.parent.mesh, missing_points, *fun_args)


                #Append newly computed point to coupling matrix, update bookkeeping
                self.points = np.vstack((self.points, missing_points))
                self.matrix = np.vstack((self.matrix, new_matrix_elems))

                #Re-compute indices of queried points, now that all should exist
                p_existing_point_idx, m_existing_point_idx = np.where((self.points == points[:, None]).all(axis=-1))

            return self.matrix[m_existing_point_idx]
