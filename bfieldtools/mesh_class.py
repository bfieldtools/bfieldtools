'''
Contains a class used for wrapping a mesh (in the form of a Trimesh object) together with
some convinient functions and properties.

'''

from time import time
import pickle

from mayavi import mlab
import trimesh
import numpy as np
from psutil import virtual_memory

from . import utils
from .laplacian_mesh import laplacian_matrix, mass_matrix
from .mutual_inductance_mesh import self_inductance_matrix


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
        Compute and return mutual inductance matrix. If mesh consists of multiple separate sub-meshes, compute these separately.
        '''

        #Available RAM in Gigabytes
        mem = virtual_memory().available >> 30

        #Estimate of memory use
        mem_per_vertex = 8 / 2000

        n_chunks = int(np.ceil(mem_per_vertex / mem * len(self.mesh.vertices)))

        print('Computing inductance matrix in %d chunks since %d GiB memory is available...'%(n_chunks, mem))

        start = time()

        inductance = self_inductance_matrix(self.mesh, Nchunks=n_chunks)

        duration = time() - start
        print('Inductance matrix computation took %.2f seconds.'%duration)

        return inductance


    @LazyProperty
    def resistance(self, resistivity=1.68*1e-8, thickness=1e-4):
        '''
        Compute and return resistance/resistivity matrix using Laplace matrix.
        Default resistivity set to that of copper.
        NB! For now, the resistivity and thickness values are set in stone.
        To continue using @LazyProperty, they could be moved into class attributes.
        Alternatively, this LazyProperty could be turned into a method.
        '''

        #Flip sign
        negative_laplacian = -1 * self.laplacian.todense()

        #Compensate for rank n-1 by adding offset
        negative_laplacian += np.ones(negative_laplacian.shape)/negative_laplacian.shape[0]

        resistance = resistivity / thickness * negative_laplacian

        #Set boundary vertices to zero
        resistance[self.boundary_verts, :][:, self.boundary_verts] = 0

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
