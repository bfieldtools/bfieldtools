from mayavi import mlab
import trimesh
import numpy as np
from psutil import virtual_memory
from time import time

from . import utils
from .laplacian_mesh import laplacian_matrix, mass_matrix
from .mutual_inductance_mesh import self_inductance_matrix, mutual_inductance_matrix


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

    def __init__(self, verts=None, tris=None, mesh_file=None):

        if mesh_file: #First, check if mesh_file passed
            self.mesh = trimesh.load(mesh_file, process=False)


        else: #Fall back on verts, tris parameters
            if isinstance(verts, type(None)) or isinstance(tris, type(None)):
                ValueError('You must provide either verts and tris or a mesh file')
            self.mesh = trimesh.Trimesh(verts, tris, process=False)

        self.verts = self.mesh.vertices
        self.tris = self.mesh.faces

        self.tri_areas = self.mesh.area_faces
        self.tri_normals = self.mesh.face_normals


        #Useful mesh metrics etc, more can be added
        self.dual_areas = utils.dual_areas(self.tris, self.tri_areas)

        self.boundary_verts, self.inner_verts,\
        self.boundary_tris, self.inner_tris = utils.find_mesh_boundaries(self.verts,
                                                                         self.tris,
                                                                         self.mesh.edges)


    @LazyProperty
    def laplacian(self):
        '''
        Compute and return surface laplacian matrix.
        '''
        laplacian = laplacian_matrix(self.verts, self.tris,
                                     self.tri_normals,
                                     self.tri_areas)

        return laplacian


    @LazyProperty
    def mass(self):
        '''
        Compute and return mesh mass matrix.
        '''
        mass = mass_matrix(self.verts, self.tris, self.tri_areas, self.dual_areas)

        return mass


    @LazyProperty
    def inductance(self):
        '''
        Compute and return mutual inductance matrix. If mesh consists of multiple separate sub-meshes, compute these separately.
        '''

        #Available RAM in Gigabytes
        mem = virtual_memory().available >> 30

        #Estimate of memory use
        mem_per_vertex = 6 / 2000

        Nchunks = int(np.ceil(mem_per_vertex / mem * len(self.verts)))

        print('Computing inductance matrix in %d chunks since %d GiB memory is available...'%(Nchunks, mem))

        start = time()
#        #If mesh corresponds of many submeshes, compute these separately to save memory
#        if self.mesh.body_count > 1:
#
#
#            inductance = np.zeros((len(self.mesh.vertices), len(self.mesh.vertices)))
#
#            #Split into separate sub-bodies
#            submeshes = self.mesh.split(only_watertight=False)
#
#            n_submeshes = len(submeshes)
#
#            #Determine how submesh vertex indices correspond to full mesh vertex indices
#            vertex_lookup = []
#
#            for i in range(n_submeshes):
#                vertex_lookup.append([])
#                for vert in submeshes[i].vertices:
#                    vertex_lookup[i].append(np.where((self.mesh.vertices == vert).all(axis=1) == True)[0][0])
#
#            #Loop through block matrix components
#            for i in range(n_submeshes):
#                for j in range(i, n_submeshes):
#                    if i==j:
#                        sub_block = self_inductance_matrix(submeshes[i].vertices, submeshes[i].faces)
#                    else:
#                        sub_block = mutual_inductance_matrix(submeshes[i].vertices, submeshes[i].faces,
#                                                 submeshes[j].vertices, submeshes[j].faces)
#
#                    #Assign to full matrix
#                    inductance[np.asarray(vertex_lookup[i])[:,None], vertex_lookup[j]] = sub_block
#
#                    #Fill in lower triangle, matrix is symmetric
#                    if i != j:
#                        inductance[np.asarray(vertex_lookup[j])[:,None], vertex_lookup[i]] = sub_block
#
#        else:

        inductance = self_inductance_matrix(self.verts, self.tris, Nchunks=Nchunks)

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

        R = resistivity / thickness
        resistance = R * self.laplacian.todense()

        #Set boundary vertices to zero
        resistance[self.boundary_verts, :][:, self.boundary_verts] = 0

        return resistance


    def plot_mesh(self):
        '''
        Simply plot the mesh surface in mayavi.
        '''

        mesh = mlab.triangular_mesh(*self.verts.T, self.tris,
                                    representation='wireframe', color=(0, 0, 0))

        return mesh
