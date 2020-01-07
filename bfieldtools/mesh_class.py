'''
Contains a class used for wrapping a mesh (in the form of a Trimesh object) together with
some convenient functions and properties.

'''

from time import time
import pickle

from mayavi import mlab
import trimesh
import numpy as np
from psutil import virtual_memory

from . import utils
from .mesh_calculus import laplacian_matrix, mass_matrix
from .mesh_inductance import self_inductance_matrix
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


class MeshWrapper:
    '''
    Class that is used for surface mesh field calculations, e.g. coil design.
    Computation functions are typically external functions that are called
    using lazy properties.

    The mesh surface can consist of a single contiguous surface or several separate
    surfaces within a single mesh object.

    '''

    def __init__(self, verts=None, tris=None, mesh_file=None, mesh_obj=None, process=False, fix_normals=True):
        '''
        Initialize MeshWrapper object.
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


        #Useful mesh metrics etc, more can be added
        self.dual_areas = utils.dual_areas(self.mesh.faces, self.mesh.area_faces)

        self.boundary_verts, self.inner_verts,\
        self.boundary_tris, self.inner_tris = utils.find_mesh_boundaries(self.mesh.vertices,
                                                                         self.mesh.faces,
                                                                         self.mesh.edges)

        self.B_coupling = CouplingMatrix(self, magnetic_field_coupling)
        self.U_coupling = CouplingMatrix(self, scalar_potential_coupling)
        self.A_coupling = CouplingMatrix(self, vector_potential_coupling)


    @LazyProperty
    def laplacian(self):
        '''
        Compute and return surface laplacian matrix.

        '''
        laplacian = laplacian_matrix(self.mesh)

        return laplacian


    @LazyProperty
    def mass(self):
        '''
        Compute and return mesh mass matrix.

        '''
        mass = mass_matrix(self.mesh, self.dual_areas)

        return mass


    @LazyProperty
    def inductance(self):
        '''
        Compute and return mutual inductance matrix.

        '''

        #Available RAM in megabytes
        mem = virtual_memory().available >> 20


        #Estimate of memory usage in megabytes for a single chunk, when quad_degree=2 (very close with quad_degree=1)
        mem_use = 0.04 * len(self.mesh.vertices)**1.75

        print('Estimating %d MiB required for %d vertices...'%(mem_use, len(self.mesh.vertices)))

        #Chunk computation so that available memory is sufficient
        n_chunks = int(np.ceil(mem_use/mem))

        print('Computing inductance matrix in %d chunks since %d MiB memory is available...'%(n_chunks, mem))

        start = time()

        inductance = self_inductance_matrix(self.mesh, Nchunks=n_chunks)

        duration = time() - start
        print('Inductance matrix computation took %.2f seconds.'%duration)

        return inductance


    def resistance(self, resistivity=1.68*1e-8, thickness=1e-4):
        '''
        Compute and return resistance/resistivity matrix using Laplace matrix.
        Default resistivity set to that of copper.
        Parameters
        ----------
        Resistivity: float
            Resistivity value in Ohm/meter
        Thickness: float
            Thickness of surface. NB! Must be small in comparison to observation distance

        Returns
        -------
        R: array
            Resistance matrix

        '''

        #Flip sign
        negative_laplacian = -1 * self.laplacian.todense()

        #Compensate for rank n-1 by adding offset
        negative_laplacian += np.ones(negative_laplacian.shape)/negative_laplacian.shape[0]

        resistance = resistivity / thickness * negative_laplacian

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
        Save the MeshWrapper object using a pickled Python file

        Parameters
        ----------
        target_file: str
            File name or file object to save to

        """

        pickle.dump(obj=self, file=open(target_file, 'wb'))


def save_pickle(obj, target_file):
    """
    Save the MeshWrapper object using a pickled Python file

    Parameters
    ----------
    obj: object to save to file
    target_file: str
        file name or file object to save to

    """

    pickle.dump(obj=obj, file=open(target_file, 'wb'), protocol=-1)


def load_pickle(target_file):
    """
    Load pickled MeshWrapper object from file

    Parameters
    -----------
    target_file: str
        File name or file object to load from
    Returns
    -------
    obj: loaded MeshWrapper object


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
