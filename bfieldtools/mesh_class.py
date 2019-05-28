from mayavi import mlab
import trimesh
import numpy as np

import utils
from laplacian_mesh import laplacian_matrix, mass_matrix
from mutual_inductance_mesh import self_inductance_matrix, mutual_inductance_matrix


class LazyProperty():
    '''
    Implementation of lazily loading properties, see
    http://blog.pythonisito.com/2008/08/lazy-descriptors.html
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


class ToBeNamed:
    '''
    Class that is used for surface mesh field calculations, e.g. coil design.
    Computation functions are typically external functions that are called
    using method wrappers.

    The mesh surface should typically consists of a single contiguous surface,
    although multi-surface meshes should also work. To combine multiple objects
    of this class, use the MultiSurface class (ToBeImplemented).
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
        Compute and return surface laplacian matrix as well as mass matrix.
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
        Compute and return mutual inductance matrix.
        '''

        #If mesh corresponds of many submeshes, compute these separately to save memory
        if self.mesh.body_count > 1:


            inductance = np.zeros((len(self.mesh.vertices), len(self.mesh.vertices)))

            #Split into separate sub-bodies
            submeshes = self.mesh.split(only_watertight=False)

            n_submeshes = len(submeshes)

            #Determine how submesh vertex indices correspond to full mesh vertex indices
            vertex_lookup = []

            for i in range(n_submeshes):
                vertex_lookup.append([])
                for vert in submeshes[i].vertices:
                    vertex_lookup[i].append(np.where((self.mesh.vertices == vert).all(axis=1) == True)[0][0])

            #Loop through block matrix components
            for i in range(n_submeshes):
                for j in range(i, n_submeshes):
                    if i==j:
                        sub_block = self_inductance_matrix(submeshes[i].vertices, submeshes[i].faces)
                    else:
                        sub_block = mutual_inductance_matrix(submeshes[i].vertices, submeshes[i].faces,
                                                 submeshes[j].vertices, submeshes[j].faces)

                    #Assign to full matrix
                    inductance[np.asarray(vertex_lookup[i])[:,None], vertex_lookup[j]] = sub_block

                    if i != j:
                        inductance[np.asarray(vertex_lookup[j])[:,None], vertex_lookup[i]] = sub_block


#            #Fill in lower triangle
#            i_lower = np.tril_indices(inductance.shape[0], -1)
#            inductance[i_lower] = inductance.T[i_lower]
        else:
            inductance = self_inductance_matrix(self.verts, self.tris)

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


    #%% Compact example of design of a biplanar coil

if __name__ == '__main__':

    import numpy as np
    from utils import cylinder_points
    from magnetic_field_mesh import compute_C
    from coil_optimize import optimize_streamfunctions
    import matplotlib.pyplot as plt


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1e2


    #Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(file_obj='./example_meshes/10x10_plane_hires.obj', process=False)

    planemesh.apply_scale(scaling_factor)

    #Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 4, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                                 planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                             planemesh.faces, process=False)

    joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = ToBeNamed(verts=joined_planes.vertices, tris=joined_planes.faces)

    #%% Set up target and stray field points

    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 0]) * scaling_factor

    sidelength = 2 * scaling_factor
    n = 6
    xx = np.linspace(-sidelength/2, sidelength/2, n)
    yy = np.linspace(-sidelength/2, sidelength/2, n)
    zz = np.linspace(-sidelength/2, sidelength/2, n)
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    target_points = np.array([x, y, z]).T

    #Turn cube into sphere by rejecting points "in the corners"
    target_points = target_points[np.linalg.norm(target_points, axis=1) < sidelength/2]  + center



#    #Here, the stray field points are on a spherical surface
    stray_radius = 20 * scaling_factor
#    stray_length = 20 * scaling_factor
#
#    stray_points = cylinder_points(radius=stray_radius,
#                                   length = stray_length,
#                                   nlength = 5,
#                                   nalpha = 30,
#                                   orientation=np.array([1, 0, 0]))
#
    stray_points_mesh = trimesh.creation.icosphere(subdivisions=2, radius=stray_radius)
    stray_points = stray_points_mesh.vertices + center

    n_stray_points = len(stray_points)



    #%% Compute C matrices that are used to compute the generated magnetic field
    coil.C = compute_C(coil.mesh, target_points)
    coil.strayC = compute_C(coil.mesh, stray_points)

    #%% Specify target field and run solver

    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function
    target_field = np.ones(target_points.shape[0], )


    # The tolerance parameter will determine the spatial detail of the coil.
    # Smaller tolerance means better but more intricate patterns. Too small values
    # will not be solveable.
    tolerance = 0.3

    I, sol = optimize_streamfunctions(coil, target_field,
                                 target_axis=0,
                                 target_error={'on_axis':0.01, 'off_axis':0.01, 'stray':0.01},
                                 laplacian_smooth=0,
                                 tolerance=tolerance)


    #%% Plot coil windings and target points

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(480, 480))
    mlab.clf()

    surface = mlab.pipeline.triangular_mesh_source(*coil.verts.T, coil.tris,scalars=I)

    windings = mlab.pipeline.contour_surface(surface, contours=8)


    B_target = np.vstack((coil.C[:, :, 0].dot(I),
                          coil.C[:, :, 1].dot(I),
                          coil.C[:, :, 2].dot(I))).T


    mlab.quiver3d(*target_points.T, *B_target.T)


    #%% Plot field falloff on two axes

    plt.figure()

    z1 = np.linspace(0, 30, 31) * scaling_factor

    x1 = y1 = np.zeros_like(z1)

    line1_points = np.vstack((x1, y1, z1)).T

    line1_C = compute_C(coil.mesh, r=line1_points)

    B_line1 = np.vstack((line1_C[:, :, 0].dot(I), line1_C[:, :, 1].dot(I), line1_C[:, :, 2].dot(I))).T

    plt.semilogy(z1 / scaling_factor, np.linalg.norm(B_line1, axis=1)/np.mean(np.abs(target_field)), label='Z')

    y2 = np.linspace(0, 30, 31) * scaling_factor

    z2 = x2 = np.zeros_like(y2)

    line2_points = np.vstack((x2, y2, z2)).T

    line2_C = compute_C(coil.mesh, r=line2_points)

    B_line2 = np.vstack((line2_C[:, :, 0].dot(I), line2_C[:, :, 1].dot(I), line2_C[:, :, 2].dot(I))).T

    plt.semilogy(y2 / scaling_factor, np.linalg.norm(B_line2, axis=1)/np.mean(np.abs(target_field)), label='Y')
    plt.ylabel('Field amplitude (target field units)')
    plt.xlabel('Distance from origin')
    plt.grid(True, which='minor', axis='y')
    plt.grid(True, which='major', axis='y', color='k')
    plt.grid(True, which='major', axis='x')

    plt.legend()
