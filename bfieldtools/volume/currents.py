"""Complete conserved face-flux spaces, including handles and cavities.

The dual-cycle construction is sparse. The local/global decomposition is a
size-limited dense reference, intended for validation and modest meshes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import qr, svd
from scipy.sparse import coo_matrix, csr_matrix, hstack
from scipy.sparse.csgraph import connected_components

from .topology import TetraTopology, _readonly


@dataclass(frozen=True, eq=False)
class ConservedFluxBasis:
    """Geometry-independent current coefficients mapped to all oriented faces.

    ``fluxes`` has shape (nfaces, n_dofs), zero boundary rows, and independent
    conserved columns. Coefficients and face fluxes have units A. Keep this
    object separately from the Cartesian CellCurrentBasis used for assembly.

    For the augmented construction, the first len(local_edges) columns are
    curls of the selected edges; the remaining n_global columns are Euclidean
    orthonormal face-flux representatives, orthogonal to all interior edge
    curls. They are not normalized to a particular physical cut. For dual
    cycles, local_edges and n_global are None (cycles mix local/global content).
    """

    topology: TetraTopology
    fluxes: csr_matrix
    local_edges: NDArray[np.int64] | None = None
    n_global: int | None = None

    @property
    def n_dofs(self) -> int:
        return self.fluxes.shape[1]


def _interior_dual(topology: TetraTopology):
    """Check boundary manifold links and return oriented dual graph edges."""
    if not isinstance(topology, TetraTopology):
        raise TypeError("topology must be a TetraTopology")
    boundary_faces = topology.faces[topology.boundary_face_mask]
    counts = np.bincount(
        topology.face_edges[topology.boundary_face_mask].ravel(),
        minlength=len(topology.edges),
    )
    if np.any(counts[topology.boundary_edge_mask] != 2):
        raise ValueError("current basis requires a closed manifold boundary")
    # Each boundary vertex must have a single circular link. Edge counts alone
    # would accept two otherwise valid solids pinched together at one vertex.
    links = {}
    for face in boundary_faces:
        for i, vertex in enumerate(face):
            links.setdefault(int(vertex), []).append(tuple(np.delete(face, i)))
    for pairs in links.values():
        pairs = np.asarray(pairs)
        _, inverse = np.unique(pairs, return_inverse=True)
        pairs = inverse.reshape(-1, 2)
        graph = coo_matrix(
            (np.ones(len(pairs)), pairs.T),
            shape=(int(pairs.max()) + 1,) * 2,
        ).tocsr()
        if connected_components(graph, directed=False)[0] != 1:
            raise ValueError("current basis requires manifold boundary vertex links")
    interior_faces = np.flatnonzero(~topology.boundary_face_mask)
    incidence = topology.divergence[:, interior_faces].tocsc()
    if np.any(np.diff(incidence.indptr) != 2):
        raise ValueError("each interior face must join two cells")
    cells = incidence.indices.reshape(-1, 2)
    signs = incidence.data.reshape(-1, 2)
    if np.any(signs.sum(axis=1) != 0):
        raise ValueError("interior faces must have opposite cell orientations")
    return interior_faces, cells, signs


def dual_cycle_basis(topology: TetraTopology) -> ConservedFluxBasis:
    """Build a complete sparse flux basis using fundamental cell-dual cycles.

    Supports disconnected manifold solids, handles, and cavities. Each chord
    carries unit positive oriented face flux; its tree path closes the cycle.
    Storage is proportional to the total cycle lengths, which can be large.
    No dense nullspace or optional meshing backend is needed.
    """
    face_ids, cells, signs = _interior_dual(topology)
    ncells = len(topology.tetrahedra)
    union = np.arange(ncells)

    def root(node):
        while union[node] != node:
            union[node] = union[union[node]]
            node = union[node]
        return node

    adjacency = [[] for _ in range(ncells)]
    chords = []
    for edge, (a, b) in enumerate(cells):
        ra, rb = root(a), root(b)
        if ra == rb:
            chords.append(edge)
        else:
            union[rb] = ra
            adjacency[a].append((b, edge, int(signs[edge, 0])))
            adjacency[b].append((a, edge, int(signs[edge, 1])))

    parent = np.full(ncells, -1, dtype=np.int64)
    parent_edge = np.full(ncells, -1, dtype=np.int64)
    parent_sign = np.zeros(ncells, dtype=np.int8)
    depth = np.zeros(ncells, dtype=np.int64)
    ncomponents = 0
    for start in range(ncells):
        if parent[start] >= 0:
            continue
        ncomponents += 1
        parent[start] = start
        stack = [start]
        while stack:
            a = stack.pop()
            for b, edge, sign_at_a in adjacency[a]:
                if parent[b] >= 0:
                    continue
                parent[b], parent_edge[b] = a, edge
                parent_sign[b] = -sign_at_a
                depth[b] = depth[a] + 1
                stack.append(b)
    if ncomponents != topology.n_components:
        raise ValueError("conductor components must be connected through faces")

    rows, columns, values = [], [], []
    for column, chord in enumerate(chords):
        rows.append(int(face_ids[chord]))
        columns.append(column)
        values.append(1)
        a, b = cells[chord]
        # Cancel the chord's balance at a by routing flux from a towards b.
        balance = -int(signs[chord, 0])
        while a != b:
            if depth[a] >= depth[b]:
                rows.append(int(face_ids[parent_edge[a]]))
                columns.append(column)
                values.append(balance * int(parent_sign[a]))
                a = parent[a]
            else:
                rows.append(int(face_ids[parent_edge[b]]))
                columns.append(column)
                values.append(-balance * int(parent_sign[b]))
                b = parent[b]
    fluxes = coo_matrix(
        (np.asarray(values, dtype=np.int64), (rows, columns)),
        shape=(len(topology.faces), len(chords)),
    ).tocsr()
    return ConservedFluxBasis(topology, fluxes)


def augmented_flux_basis(
    topology: TetraTopology, *, rtol: float = 1e-10, max_dense_entries: int = 4_000_000
) -> ConservedFluxBasis:
    """Independent interior-edge curls plus all missing global current modes.

    Pivoted QR removes all curl-free edge potentials, including cavity gauge
    modes. A projected dual-cycle basis supplies the handle complement. The
    result spans exactly the same space as dual_cycle_basis(). This dense
    reference refuses large intermediates; use dual_cycle_basis for sparse
    construction. rtol is relative to the largest QR pivot for local rank and
    to unit scale for the projected orthonormal conserved space.
    """
    if not np.isfinite(rtol) or not 0 < rtol < 1:
        raise ValueError("rtol must be finite and between zero and one")
    if (
        isinstance(max_dense_entries, (bool, np.bool_))
        or not isinstance(max_dense_entries, (int, np.integer))
        or max_dense_entries < 1
    ):
        raise ValueError("max_dense_entries must be a positive integer")
    # Check allocation sizes before constructing the possibly long dual cycles.
    face_ids, _, _ = _interior_dual(topology)
    edge_ids = np.flatnonzero(~topology.boundary_edge_mask)
    ncycles = len(face_ids) - len(topology.tetrahedra) + topology.n_components
    if max(len(face_ids) * len(edge_ids), len(face_ids) * ncycles) > max_dense_entries:
        raise ValueError(
            "dense reference exceeds max_dense_entries; use dual_cycle_basis"
        )
    cycles = dual_cycle_basis(topology)
    local = topology.curl[face_ids][:, edge_ids].toarray().astype(float)
    q, r, pivots = qr(local, mode="economic", pivoting=True)
    diagonal = np.abs(np.diag(r))
    rank = (
        int(np.count_nonzero(diagonal > rtol * diagonal.max())) if diagonal.size else 0
    )
    local_edges = edge_ids[pivots[:rank]]
    _readonly(local_edges)
    n_global = ncycles - rank
    if n_global < 0:
        raise ValueError("inconsistent numerical ranks; reduce rtol")
    z = cycles.fluxes[face_ids].toarray().astype(float)
    z, _ = qr(z, mode="economic")
    complement = z - q[:, :rank] @ (q[:, :rank].T @ z)
    u, singular, _ = svd(complement, full_matrices=False)
    if np.count_nonzero(singular > rtol) != n_global:
        raise ValueError("inconsistent numerical ranks; adjust rtol")
    global_fluxes = coo_matrix(
        (
            u[:, :n_global].ravel(),
            (
                np.repeat(face_ids, n_global),
                np.tile(np.arange(n_global), len(face_ids)),
            ),
        ),
        shape=(len(topology.faces), n_global),
    ).tocsr()
    fluxes = hstack((topology.curl[:, local_edges], global_fluxes), format="csr")
    return ConservedFluxBasis(topology, fluxes, local_edges, n_global)
