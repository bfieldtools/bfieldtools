"""Complete conserved face-flux spaces, including handles and cavities.

Global modes are the right nullspace of the sparse stacked divergence and
transpose-curl operators. QR/SVD use a size-limited dense reference solve.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import qr, svd
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack
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
    curls. They are not normalized to a particular physical cut.
    """

    topology: TetraTopology
    fluxes: csr_matrix
    local_edges: NDArray[np.int64] | None = None
    n_global: int | None = None

    @property
    def n_dofs(self) -> int:
        return self.fluxes.shape[1]


def _interior_faces(topology: TetraTopology):
    """Check manifold boundaries and face connectivity; return interior face IDs."""
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
    graph = coo_matrix(
        (np.ones(len(cells)), cells.T),
        shape=(len(topology.tetrahedra),) * 2,
    ).tocsr()
    if connected_components(graph, directed=False)[0] != topology.n_components:
        raise ValueError("conductor components must be connected through faces")
    return interior_faces


def augmented_flux_basis(
    topology: TetraTopology, *, rtol: float = 1e-10, max_dense_entries: int = 4_000_000
) -> ConservedFluxBasis:
    """Independent interior-edge curls plus stacked-nullspace global modes.

    With boundary fluxes fixed to zero, solve D_i H = 0 and C_i.T H = 0.
    These are Euclidean face-flux harmonic representatives (identity metric),
    not geometry-weighted RT0 harmonic fields. Pivoted QR selects independent
    local curls, removing all potential gauges including cavity redundancies.
    The combined columns span the complete conserved face-flux space.

    The stack is sparse, but this reference uses dense QR/SVD. max_dense_entries
    bounds each dense input and the square right-singular-vector workspace;
    it is not a total-memory bound. rtol is relative to the largest QR pivot
    for local rank and the largest singular value for stacked rank.
    """
    if not np.isfinite(rtol) or not 0 < rtol < 1:
        raise ValueError("rtol must be finite and between zero and one")
    if (
        isinstance(max_dense_entries, (bool, np.bool_))
        or not isinstance(max_dense_entries, (int, np.integer))
        or max_dense_entries < 1
    ):
        raise ValueError("max_dense_entries must be a positive integer")
    face_ids = _interior_faces(topology)
    edge_ids = np.flatnonzero(~topology.boundary_edge_mask)
    divergence = topology.divergence[:, face_ids].astype(float)
    local = topology.curl[face_ids][:, edge_ids].astype(float)
    stacked = vstack((divergence, local.T), format="csr")
    nfaces = len(face_ids)
    if (
        max(nfaces * len(edge_ids), stacked.shape[0] * nfaces, nfaces**2)
        > max_dense_entries
    ):
        raise ValueError("dense reference exceeds max_dense_entries")
    _, r, pivots = qr(local.toarray(), mode="economic", pivoting=True)
    diagonal = np.abs(np.diag(r))
    rank = (
        int(np.count_nonzero(diagonal > rtol * diagonal.max())) if diagonal.size else 0
    )
    local_edges = edge_ids[pivots[:rank]]
    _readonly(local_edges)
    dimension = nfaces - len(topology.tetrahedra) + topology.n_components
    expected_global = dimension - rank
    # Economic SVD retains the entire right space for tall stacks. Wide
    # stacks need full Vh to include the implicit zero singular values.
    _, singular, vh = svd(stacked.toarray(), full_matrices=stacked.shape[0] < nfaces)
    stacked_rank = (
        int(np.count_nonzero(singular > rtol * singular.max())) if singular.size else 0
    )
    harmonic = vh[stacked_rank:].T
    n_global = harmonic.shape[1]
    if expected_global < 0 or n_global != expected_global:
        raise ValueError("inconsistent numerical ranks; adjust rtol")
    global_fluxes = coo_matrix(
        (
            harmonic.ravel(),
            (
                np.repeat(face_ids, n_global),
                np.tile(np.arange(n_global), len(face_ids)),
            ),
        ),
        shape=(len(topology.faces), n_global),
    ).tocsr()
    fluxes = hstack((topology.curl[:, local_edges], global_fluxes), format="csr")
    return ConservedFluxBasis(topology, fluxes, local_edges, n_global)
