"""Oriented mesh connectivity and a boundary-compatible edge-potential gauge.

Geometry lives in TetraMesh. This module contains integer connectivity and
sparse incidence operators; it does not assemble material or magnetic matrices.
"""

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree

from .mesh import TETRA_EDGES, TETRA_FACES, TetraMesh


def _readonly(*arrays):
    for array in arrays:
        array.flags.writeable = False


@dataclass(frozen=True, eq=False)
class TetraTopology:
    """Connectivity and oriented incidence operators for a tetrahedral mesh.

    Construct with from_mesh() or build_topology(). Frozen fields prevent
    reassignment; SciPy matrix storage remains mutable.

    Global edges point from the smaller vertex index to the larger. Global
    faces have sorted vertex indices, with orientation (a, b, c). Local edges
    follow TETRA_EDGES; local faces follow the outward TETRA_FACES ordering.

    ``gradient`` maps nodal values to edge differences, ``curl`` maps edge
    integrals to face fluxes, and ``divergence`` maps face fluxes to cell
    balances. These are incidence matrices, without metric factors. In
    particular, divergence gives integrated balance, not divergence density.

    Unused mesh vertices are allowed and have component label -1. Face, edge,
    and vertex boundary masks refer to the boundary of the volume mesh.
    """

    mesh: TetraMesh
    edges: np.ndarray
    cell_edges: np.ndarray
    cell_edge_signs: np.ndarray
    faces: np.ndarray
    cell_faces: np.ndarray
    cell_face_signs: np.ndarray
    face_edges: np.ndarray
    boundary_face_mask: np.ndarray
    boundary_edge_mask: np.ndarray
    boundary_vertex_mask: np.ndarray
    active_vertex_mask: np.ndarray
    vertex_components: np.ndarray
    cell_components: np.ndarray
    n_components: int

    @classmethod
    def from_mesh(cls, mesh):
        """Construct through the separate build_topology() algorithm."""
        return build_topology(mesh)

    @cached_property
    def gradient(self):
        """Cached edge–vertex incidence; see gradient_matrix()."""
        return gradient_matrix(self)

    @cached_property
    def curl(self):
        """Cached face–edge incidence; see curl_matrix()."""
        return curl_matrix(self)

    @cached_property
    def divergence(self):
        """Cached cell–face incidence; see divergence_matrix()."""
        return divergence_matrix(self)


@dataclass(frozen=True, eq=False)
class TreeCotreeGauge:
    """Independent edge-potential representation after boundary and tree gauging.

    tree_edges and retained_edges index the topology's global edge array.
    Boundary edges belong to neither set. prolongation maps independent
    coefficients to all edge coefficients; curl maps them to all face fluxes.
    n_dofs is the number of retained edges. Construct with from_topology() or
    tree_cotree_gauge().
    Frozen fields prevent reassignment; SciPy matrix storage remains mutable.
    """

    tree_edges: np.ndarray
    retained_edges: np.ndarray
    n_dofs: int
    prolongation: csr_matrix
    curl: csr_matrix

    @classmethod
    def from_topology(cls, topology):
        """Construct through the separate tree_cotree_gauge() algorithm."""
        return tree_cotree_gauge(topology)


def build_topology(mesh):
    """Build oriented connectivity and sparse incidence matrices once.

    This function does not choose a gauge or assume a solid-ball topology.
    """
    if not isinstance(mesh, TetraMesh):
        raise TypeError("mesh must be a TetraMesh")
    local_edges = mesh.tetrahedra[:, TETRA_EDGES]
    edges, inverse = np.unique(
        np.sort(local_edges.reshape(-1, 2), axis=1),
        axis=0,
        return_inverse=True,
    )
    cell_edges = inverse.reshape(-1, 6)
    cell_edge_signs = np.where(
        local_edges[:, :, 0] < local_edges[:, :, 1], 1, -1
    ).astype(np.int8)

    local_faces = mesh.tetrahedra[:, TETRA_FACES].reshape(-1, 3)
    faces, inverse, counts = np.unique(
        np.sort(local_faces, axis=1),
        axis=0,
        return_inverse=True,
        return_counts=True,
    )
    cell_faces = inverse.reshape(-1, 4)
    # Parity of each outward face relative to its sorted global orientation.
    inversions = (
        (local_faces[:, 0] > local_faces[:, 1]).astype(np.int8)
        + (local_faces[:, 0] > local_faces[:, 2])
        + (local_faces[:, 1] > local_faces[:, 2])
    )
    cell_face_signs = (1 - 2 * (inversions % 2)).reshape(-1, 4)
    boundary_face_mask = counts == 1

    # Each sorted face has boundary (a,b) + (b,c) - (a,c).
    face_pairs = faces[:, ((0, 1), (1, 2), (0, 2))]
    # A sparse lookup avoids packing pairs into potentially overflowing keys.
    lookup = csr_matrix(
        (np.arange(len(edges)) + 1, edges.T),
        shape=(len(mesh.vertices), len(mesh.vertices)),
    )
    face_edges = (
        np.asarray(
            lookup[face_pairs[:, :, 0].ravel(), face_pairs[:, :, 1].ravel()]
        ).reshape(-1, 3)
        - 1
    ).astype(np.int64)
    boundary_edge_mask = np.zeros(len(edges), dtype=bool)
    boundary_edge_mask[face_edges[boundary_face_mask]] = True
    boundary_vertex_mask = np.zeros(len(mesh.vertices), dtype=bool)
    boundary_vertex_mask[faces[boundary_face_mask]] = True
    active_vertex_mask = np.zeros(len(mesh.vertices), dtype=bool)
    active_vertex_mask[edges] = True

    adjacency = coo_matrix(
        (np.ones(len(edges), dtype=np.int8), edges.T),
        shape=(len(mesh.vertices), len(mesh.vertices)),
    ).tocsr()
    _, labels = connected_components(adjacency, directed=False)
    active_labels, inverse = np.unique(labels[active_vertex_mask], return_inverse=True)
    n_components = len(active_labels)
    vertex_components = np.full(len(mesh.vertices), -1, dtype=np.int64)
    vertex_components[active_vertex_mask] = inverse
    cell_components = vertex_components[mesh.tetrahedra[:, 0]]
    _readonly(
        edges,
        cell_edges,
        cell_edge_signs,
        faces,
        cell_faces,
        cell_face_signs,
        face_edges,
        boundary_face_mask,
        boundary_edge_mask,
        boundary_vertex_mask,
        active_vertex_mask,
        vertex_components,
        cell_components,
    )

    return TetraTopology(
        mesh=mesh,
        edges=edges,
        cell_edges=cell_edges,
        cell_edge_signs=cell_edge_signs,
        faces=faces,
        cell_faces=cell_faces,
        cell_face_signs=cell_face_signs,
        face_edges=face_edges,
        boundary_face_mask=boundary_face_mask,
        boundary_edge_mask=boundary_edge_mask,
        boundary_vertex_mask=boundary_vertex_mask,
        active_vertex_mask=active_vertex_mask,
        vertex_components=vertex_components,
        cell_components=cell_components,
        n_components=n_components,
    )


def gradient_matrix(topology):
    """Edge–vertex incidence G: edge differences = G @ nodal values.

    Shape (nedges, nvertices); each row is -1 at the tail, +1 at the head.
    """
    return coo_matrix(
        (
            np.tile([-1, 1], len(topology.edges)).astype(np.int8),
            (np.repeat(np.arange(len(topology.edges)), 2), topology.edges.ravel()),
        ),
        shape=(len(topology.edges), len(topology.mesh.vertices)),
    ).tocsr()


def curl_matrix(topology):
    """Face–edge incidence C: oriented face fluxes = C @ edge integrals.

    Shape (nfaces, nedges); each sorted face (a,b,c) has boundary ab + bc - ac.
    """
    return coo_matrix(
        (
            np.tile([1, 1, -1], len(topology.faces)).astype(np.int8),
            (np.repeat(np.arange(len(topology.faces)), 3), topology.face_edges.ravel()),
        ),
        shape=(len(topology.faces), len(topology.edges)),
    ).tocsr()


def divergence_matrix(topology):
    """Cell–face incidence D: outward cell balances = D @ face fluxes.

    Shape (ncells, nfaces). Divide balances by cell volumes to obtain divergence
    density. The exact incidence identities are C @ G = 0 and D @ C = 0.
    """
    return coo_matrix(
        (
            topology.cell_face_signs.ravel(),
            (
                np.repeat(np.arange(len(topology.mesh.tetrahedra)), 4),
                topology.cell_faces.ravel(),
            ),
        ),
        shape=(len(topology.mesh.tetrahedra), len(topology.faces)),
    ).tocsr()


def tree_cotree_gauge(topology):
    """Choose independent edge curls for solid-ball conductor components.

    Set boundary-edge coefficients to zero, contract each component's boundary
    to one root, and build a spanning forest joining its interior vertices to
    that root. All remaining nonboundary edges are the independent cotree
    unknowns. Parallel edges to a root and root self-loops are supported.

    Handles, cavities, and nonmanifold boundaries need additional treatment
    and raise ValueError. Geometry must satisfy TetraMesh's conformity contract.
    """
    if not isinstance(topology, TetraTopology):
        raise TypeError("topology must be a TetraTopology")
    _check_supported_boundary(topology)
    interior = topology.active_vertex_mask & ~topology.boundary_vertex_mask
    ninterior = np.count_nonzero(interior)
    quotient = np.full(len(topology.mesh.vertices), -1, dtype=np.int64)
    quotient[interior] = np.arange(ninterior)
    boundary = topology.boundary_vertex_mask
    quotient[boundary] = ninterior + topology.vertex_components[boundary]
    candidates = np.flatnonzero(~topology.boundary_edge_mask)
    pairs = np.sort(quotient[topology.edges[candidates]], axis=1)
    nonloops = pairs[:, 0] != pairs[:, 1]
    # Multiple edges can join one interior vertex to different vertices of
    # the contracted boundary. Choose the first global edge for that pair.
    pairs, first = np.unique(pairs[nonloops], axis=0, return_index=True)
    edge_ids = candidates[nonloops][first]
    nquotient = ninterior + topology.n_components
    graph = coo_matrix(
        (
            np.concatenate((edge_ids + 1, edge_ids + 1)),
            (
                np.concatenate((pairs[:, 0], pairs[:, 1])),
                np.concatenate((pairs[:, 1], pairs[:, 0])),
            ),
        ),
        shape=(nquotient, nquotient),
    ).tocsr()
    # Positive edge IDs are deterministic weights and also identify the
    # selected mesh edges. The choice of spanning tree has no physical role.
    tree_edges = np.sort(minimum_spanning_tree(graph).data.astype(np.int64) - 1)
    if len(tree_edges) != ninterior:
        raise ValueError("interior vertices do not connect to the boundary roots")
    retained = ~topology.boundary_edge_mask.copy()
    retained[tree_edges] = False
    retained_edges = np.flatnonzero(retained)
    n_dofs = len(retained_edges)

    # Compare with the complete dual-graph conserved-flux space. In this
    # supported topology the relative curl has precisely this dimension.
    interior_faces = ~topology.boundary_face_mask
    incidence = topology.divergence[:, interior_faces].astype(np.int64)
    ndual, _ = connected_components(incidence @ incidence.T, directed=False)
    expected = np.count_nonzero(interior_faces) - len(topology.mesh.tetrahedra) + ndual
    if ndual != topology.n_components or n_dofs != expected:
        raise ValueError(
            "mesh topology requires additional global patterns or constraints; "
            "the solid-ball tree gauge is insufficient"
        )
    _readonly(tree_edges, retained_edges)

    prolongation = coo_matrix(
        (np.ones(n_dofs, dtype=np.int8), (retained_edges, np.arange(n_dofs))),
        shape=(len(topology.edges), n_dofs),
    ).tocsr()
    return TreeCotreeGauge(
        tree_edges=tree_edges,
        retained_edges=retained_edges,
        n_dofs=n_dofs,
        prolongation=prolongation,
        curl=topology.curl[:, retained_edges].tocsr(),
    )


def _check_supported_boundary(topology):
    faces = topology.faces[topology.boundary_face_mask]
    face_edges = topology.face_edges[topology.boundary_face_mask]
    counts = np.bincount(face_edges.ravel(), minlength=len(topology.edges))
    if np.any(counts[topology.boundary_edge_mask] != 2):
        raise ValueError("tree gauge requires a closed manifold boundary")
    boundary_vertices = np.flatnonzero(topology.boundary_vertex_mask)
    boundary_edges = topology.edges[topology.boundary_edge_mask]
    adjacency = coo_matrix(
        (np.ones(len(boundary_edges)), boundary_edges.T),
        shape=(len(topology.mesh.vertices), len(topology.mesh.vertices)),
    ).tocsr()[boundary_vertices][:, boundary_vertices]
    _, labels = connected_components(adjacency, directed=False)
    components = topology.vertex_components[boundary_vertices]
    for component in range(topology.n_components):
        if len(np.unique(labels[components == component])) != 1:
            raise ValueError(
                "tree gauge requires one connected boundary per conductor; "
                "cavities need additional gauge constraints"
            )
        nv = np.count_nonzero(components == component)
        ne = np.count_nonzero(
            topology.vertex_components[boundary_edges[:, 0]] == component
        )
        nf = np.count_nonzero(topology.vertex_components[faces[:, 0]] == component)
        if nv - ne + nf != 2:
            raise ValueError(
                "tree gauge currently supports genus-zero boundaries; "
                "handles need additional global current patterns"
            )
