"""
This module contains miscellaneous utility functions used across bfieldtools.
"""

__all__ = [
    "MeshProjection",
    "combine_meshes",
    "cylinder_points",
    "dual_areas",
    "find_mesh_boundaries",
    "fix_normals",
    "get_line_quad_points",
    "get_quad_points",
    "inner2vert",
    "load_example_mesh",
    "tri_normals_and_areas",
    "vert2inner",
]


import os
import numpy as np
import quadpy
import pkg_resources
import trimesh


def combine_meshes(meshes):
    """
    Combine two or more non-overlapping Trimesh meshes without any dependency
    requirements. For more demanding applications, use Trimesh boolean operations

    Parameters
    ----------
    meshes: list or tuple
        Each element should be a Trimesh mesh

    Returns
    -------
    combined_mesh: Trimesh mesh
    """

    N_meshes = len(meshes)

    vertices = np.zeros((0, 3))
    faces = np.zeros((0, 3))

    for idx in range(N_meshes):

        faces = np.vstack((faces, meshes[idx].faces + len(vertices)))
        vertices = np.vstack((vertices, meshes[idx].vertices))

    combined_mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    return combined_mesh


def tri_normals_and_areas(r, tri):
    """ Get triangle normals and areas from vertices (r) and
        triangle indices (tri)
    """
    n = np.cross(
        r[tri[:, 1], :] - r[tri[:, 0], :], r[tri[:, 2], :] - r[tri[:, 0], :], axis=1
    )
    a = np.linalg.norm(n, axis=1)
    n = n / a[:, None]
    a /= 2
    return n, a


def get_quad_points(
    verts, tris, method="seven_point", index=None, return_ref_coords=False
):
    """ Get quad points and weights from quadrature rules implemented in
        quadpy

        Parameters
        ----------
        verts: array-like [Nverts x 3]
        tris: array-like [Ntris x 3]

        Returns
        -------
        w: array-like  (Nquad, )
            quadrature weights
        qp: array-like (Ntris, Nquad)
            quadrature points in each triangle

    """
    # methods = [k for k in quadpy.triangle.__dict__.keys()]# if k[0].isupper()]
    methods = [k for k in quadpy.triangle.__all__]
    if method in methods:
        try:
            rule = quadpy.triangle.__dict__[method]()
        except (TypeError) as error:
            if index is not None:
                rule = quadpy.triangle.__dict__[method](index)
            else:
                print("The method requires index (check quadpy documentation)")
                raise error
    else:
        # raise ValueError('method: '+method+' not in the available list of methods: ' + methods)
        raise ValueError(
            "method: "
            + method
            + " not in the available list of quadrature methods (check quadpy.triangle.__all__)"
        )

    x = rule.points[:, 0:2]
    w = rule.weights

    qp = np.zeros((tris.shape[0], len(w), 3))
    for i, t in enumerate(tris):
        p0 = verts[t[0]]
        p1 = verts[t[1]]
        p2 = verts[t[2]]
        B = np.array([p1 - p0, p2 - p0])

        qp[i] = x @ B + p0

    if return_ref_coords:
        return w, qp, x
    else:
        return w, qp


def get_line_quad_points(line_vertices, method="midpoint", index=None):
    """ Get quad points and weights from quadrature rules implemented in
        quadpy

        Parameters
        ----------
        line_vertices: array-like [Nverts x 3]
            Assumes vertices are connected according to index, last index connects to first

        Returns
        -------
        w: array-like  (Nquad, )
            quadrature weights
        qp: array-like (Nverts, Nquad)
            quadrature points for each edge connecting the vertices

    """
    methods = [k for k in quadpy.line_segment.__dict__.keys()]  # if k[0].isupper()]
    if method in methods:
        try:
            rule = quadpy.line_segment.__dict__[method]()
        except (TypeError) as error:
            if index is not None:
                rule = quadpy.line_segment.__dict__[method](index)
            else:
                print("The method requires index (check quadpy documentation)")
                raise error
    else:
        raise ValueError(
            "method: "
            + method
            + " not in the available list of methods (check quadpy.line_segment.__all__)"
        )

    x = rule.points
    w = rule.weights

    #    tris -> edges

    qp = np.zeros((len(line_vertices), len(w), 3))

    for i in range(len(line_vertices)):
        p0 = line_vertices[i]

        if i == len(line_vertices) - 1:
            p1 = line_vertices[0]
        else:
            p1 = line_vertices[i + 1]

        B = np.array([p1 - p0])

        qp[i] = x[:, None] @ B / 2 + B / 2 + p0

    return w, qp


def dual_areas(tris, ta):
    """ Calculate (dual) areas for each node in inds

        Dual area == area summed over the neighbouring triangles divided by 3
    """
    areas = np.zeros(np.max(tris) + 1)
    for i in range(tris.shape[0]):
        for j in range(tris.shape[1]):
            areas[tris[i, j]] += ta[i]

    return areas / 3


def find_mesh_boundaries(mesh):
    """
    Finds the open boundaries of a mesh by finding the edges that only
    belong to a single triangle. Returns an index array of inner vertices
    and triangles that do not touch the outer boundary.
    Takes edge parameter for convenience.

    Parameters
    ----------
    mesh: trimesh mesh object

    Returns
    -------
    boundaries: list of array-like

    """
    inner_vertices = np.arange(0, len(mesh.vertices))

    outline = mesh.outline(process=False)

    boundaries = []
    for i in outline.entities:
        boundaries.append(np.unique(i.points))

        inner_vertices = np.setdiff1d(inner_vertices, i.points)

    return boundaries, inner_vertices


def inner2vert(mesh, inner_vertices, holes):
    """ Linear mapping of the inner (free) weights in the stream function
        discretization to weights in all vertices

        Parameters:
            mesh: Trimesh object
            inner_vertices: list of indices of the inner vertices of the mesh
            holes: list of indices for holes in the mesh

        Returns:
            NxM sparse array, where N==mesh.vertices.shape[0]
            and M == len(inner_vertices) + len(holes)
    """
    from scipy.sparse import csr_matrix

    N = mesh.vertices.shape[0]
    M = len(inner_vertices) + len(holes)
    ii = list(inner_vertices)  # indices of inner vertices
    jj = list(np.arange(len(inner_vertices)))  # indices of inner vertices in dof
    # Hole values maps to value for each hole vertex
    for n, h in enumerate(holes):
        ii.extend(list(h))  # indices of hole indices
        jj.extend(
            [len(inner_vertices) + n] * len(h)
        )  # len(h) times index of hole in dof
    d2v = csr_matrix((np.ones(len(jj)), (ii, jj)), shape=(N, M), dtype=float)

    return d2v


def vert2inner(mesh, inner_vertices, holes):
    """ Linear mapping of the all weights in the stream function
        discretization to inner (free) weights

        Parameters:
            mesh: Trimesh object
            inner_vertices: list of indices of the inner vertices of the mesh
            holes: list of indices for holes in the mesh

        Returns:
            MxN sparse array, where N==mesh.vertices.shape[0]
            and M == len(inner_vertices) + len(holes)
    """
    from scipy.sparse import csr_matrix

    N = mesh.vertices.shape[0]
    M = len(inner_vertices) + len(holes)
    ii = list(inner_vertices)  # indices of inner vertices
    jj = list(np.arange(len(inner_vertices)))  # indices of inner vertices in dof
    vals = list(np.ones(len(inner_vertices)))
    for n, h in enumerate(holes):
        ii.extend(list(h))  # indices of hole indices
        jj.extend(
            [len(inner_vertices) + n] * len(h)
        )  # len(h) times index of hole in free values
        # Values at holes map to their average (ok, when constant boundary condition satisfied)
        vals.extend(list(np.ones(len(h)) / len(h)))
    v2d = csr_matrix((vals, (jj, ii)), shape=(M, N), dtype=float)

    return v2d


def cylinder_points(
    radius=1,
    length=1,
    nlength=10,
    alpha=360,
    nalpha=10,
    center=np.array([0, 0, 0]),
    orientation=np.array([1, 0, 0]),
):
    """
    Generate and return a set of points on a cylindrical surface.
    """
    # Create the length array
    I = np.linspace(0, length, nlength)

    # Create alpha array avoid duplication of endpoints
    # Conditional should be changed to meet your requirements
    if int(alpha) == 360:
        A = np.linspace(0, alpha, num=nalpha, endpoint=False) / 180 * np.pi
    else:
        A = np.linspace(0, alpha, num=nalpha) / 180 * np.pi

    # Calculate X and Y
    X = radius * np.cos(A)
    Y = radius * np.sin(A)

    # Tile/repeat indices so all unique pairs are present
    pz = np.tile(I, nalpha)
    px = np.repeat(X, nlength)
    py = np.repeat(Y, nlength)

    points = np.vstack((pz, px, py)).T

    # Shift to center
    shift = np.array(center) - np.mean(points, axis=0)
    points += shift

    # Orient tube to new vector

    def rotation_matrix(axis, theta):
        a = np.cos(theta / 2)
        b, c, d = -axis * np.sin(theta / 2)
        return np.array(
            [
                [
                    a * a + b * b - c * c - d * d,
                    2 * (b * c - a * d),
                    2 * (b * d + a * c),
                ],
                [
                    2 * (b * c + a * d),
                    a * a + c * c - b * b - d * d,
                    2 * (c * d - a * b),
                ],
                [
                    2 * (b * d - a * c),
                    2 * (c * d + a * b),
                    a * a + d * d - b * b - c * c,
                ],
            ]
        )

    ovec = orientation / np.linalg.norm(orientation)
    cylvec = np.array([1, 0, 0])

    if np.allclose(cylvec, ovec):
        return points

    # Get orthogonal axis and rotation
    oaxis = np.cross(ovec, cylvec)
    rot = np.arccos(np.dot(ovec, cylvec))

    R = rotation_matrix(oaxis, rot)
    return points.dot(R)


def fix_normals(mesh, origin=np.array([0, 0, 0])):
    """
    Attempts to fix face windings and normals such that normals are always "pointing out"
    from the origin.

    Parameters
    ----------
    mesh: Trimesh mesh object
    origin: array-like (3, )
        Specified from where the normals should "point out"

    Returns
    -------
    mesh: modified Trimesh object


    """

    # Dot product of all face normals and the corresponding triangle centers
    dotprods = np.sum(mesh.face_normals * (mesh.triangles_center - origin), axis=-1)
    # Flip windings for normals pointing inwards
    mesh.faces[dotprods < 0] = mesh.faces[dotprods < 0, ::-1]
    # old_cache = {'key': mesh.cache.value, ...} # Could save some old values...
    # Clear cached values
    mesh._cache.clear()
    # Could update the cache with values that don't change when flipping triangles
    # self._cache.update(old_cache)
    return mesh


def load_example_mesh(mesh_name, process=True, **kwargs):
    """
    Convenience function used load example meshes included with the package

    Parameters
    ----------
    mesh_name: string
        name of mesh, i.e. filename without extension
    process: Boolean
        Whether trimesh should process the mesh on loading
    kwargs
        Passed to trimesh object creation

    Returns
    -------
    Trimesh object
    """
    existing_files = pkg_resources.resource_listdir("bfieldtools", "example_meshes")

    # Filter according to file extension
    existing_files = [
        file
        for file in existing_files
        if file.lower().endswith(tuple(trimesh.exchange.load.mesh_formats()))
    ]

    # Remove file extension to get name
    existing_names = [os.path.splitext(file)[0] for file in existing_files]

    # Check if name exists
    if mesh_name not in existing_names:
        raise ValueError(
            "Mesh with name %s not found in example_meshes folder" % mesh_name
        )

    filename = existing_files[existing_names.index(mesh_name)]

    return trimesh.load(
        pkg_resources.resource_filename("bfieldtools", "example_meshes/" + filename),
        process=process,
        **kwargs
    )


class MeshProjection:
    """ Class for constructing projection of an arbitrary function
        to the hat functions of the given mesh
        
        The constructor initializes quadrature points, point weights
        and vertex mappings 
        
        After the construction the method "hatfunc_innerproducts" can
        be called to obtain inner products of the hat functions and 
        a given function (func)
    """

    def __init__(self, mesh, quad_degree):
        self.mesh = mesh

        weights, self.quadpoints, ref_coords = get_quad_points(
            mesh.vertices, mesh.faces, "dunavant_0" + str(quad_degree), None, True
        )
        # Mapping from coordinates of a reference triangle to barycentric coordinates
        self.M = np.array([[-1.0, -1.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        # Solve barycentric coordinates for the quad points
        # These are the shape function values for the quad points
        self.bary_coords = (
            np.hstack((ref_coords, np.ones((ref_coords.shape[0], 1)))) @ self.M.T
        )
        # Weight the shape function values by quadrature weights
        self.weights = self.bary_coords * weights[:, None]
        from scipy.sparse import csc_matrix

        Nf = len(mesh.faces)
        Nv = len(mesh.vertices)
        self.M0 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 0])), (Nf, Nv))
        self.M1 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 1])), (Nf, Nv))
        self.M2 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 2])), (Nf, Nv))

    def hatfunc_innerproducts(self, func):
        ff = np.zeros((self.quadpoints.shape[0], 3))
        for i in range(self.quadpoints.shape[1]):
            f = func(self.quadpoints[:, i, :])
            ff += self.weights[i] * f[:, None]

        ff *= self.mesh.area_faces[:, None]
        fv = ff[:, 0] @ self.M0 + ff[:, 1] @ self.M1 + ff[:, 2] @ self.M2

        return fv
