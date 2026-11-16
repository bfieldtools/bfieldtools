"""Cell-current representations and construction from oriented volume edges."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import coo_matrix, csr_matrix, isspmatrix_csr

from .mesh import TETRA_EDGES, TetraMesh

if TYPE_CHECKING:
    from .topology import TetraTopology, TreeCotreeGauge


@dataclass(frozen=True, eq=False)
class CellCurrentBasis:
    """A piecewise-constant current space on an ordered tetrahedral mesh.

    Each real CSR operator maps coefficients to one Cartesian current-density
    component in every cell; all have shape (ncells, n_dofs). Rows follow
    mesh.tetrahedra exactly. With coefficients in A and coordinates in metres,
    operators have units 1/m².

    The container holds geometry and operators without topology or gauge
    metadata. It checks layout and real finite values. Conservation, boundary
    conditions, and independence are responsibilities of the construction
    algorithm; build_edge_current_basis supplies these for the supported gauge.
    Frozen fields prevent reassignment; sparse storage remains mutable.
    """

    mesh: TetraMesh
    operators: tuple[csr_matrix, csr_matrix, csr_matrix]

    def __post_init__(self) -> None:
        if len(self.operators) != 3 or any(
            not isspmatrix_csr(operator) for operator in self.operators
        ):
            raise ValueError("operators must contain three CSR matrices")
        shape = self.operators[0].shape
        if shape[0] != len(self.mesh.tetrahedra) or any(
            operator.shape != shape for operator in self.operators
        ):
            raise ValueError("operators must share shape (ncells, n_dofs)")
        if any(
            np.iscomplexobj(operator.data) or not np.all(np.isfinite(operator.data))
            for operator in self.operators
        ):
            raise ValueError("current operators must be real and finite")

    @property
    def n_dofs(self) -> int:
        """Number of coefficient columns, independent of their construction."""
        return self.operators[0].shape[1]


def build_edge_current_basis(
    mesh: TetraMesh, topology: TetraTopology, gauge: TreeCotreeGauge
) -> CellCurrentBasis:
    """Build cell-current operators from an explicitly chosen edge gauge.

    Mesh and topology must have identical ordered connectivity. The result does
    not retain topology or gauge; keep those separately for flux interpretation.
    """
    return CellCurrentBasis(
        mesh, cell_current_operators(mesh, topology, gauge.retained_edges)
    )


def _check_connectivity(mesh: TetraMesh, topology: TetraTopology) -> None:
    if len(mesh.vertices) != topology.n_vertices or not np.array_equal(
        mesh.tetrahedra, topology.tetrahedra
    ):
        raise ValueError("mesh and topology must have identical ordered connectivity")


def cell_edge_current_vectors(
    mesh: TetraMesh, topology: TetraTopology
) -> NDArray[np.float64]:
    """Cartesian edge-current vectors in each cell, shape (ncells, 6, 3).

    Entry [T,n] is the curl of the GLOBAL function for topology.cell_edges[T,n].
    With local ordering (01,02,03,23,31,12), it is the opposing directed edge
    divided by 3V, multiplied by the potential edge's local-to-global sign.
    """
    _check_connectivity(mesh, topology)
    return (
        mesh.cell_edge_vectors[:, (np.arange(6) + 3) % 6]
        * topology.cell_edge_signs[:, :, None]
        / (3 * mesh.volumes[:, None, None])
    )


def cell_current_operators(
    mesh: TetraMesh, topology: TetraTopology, retained_edges: ArrayLike
) -> tuple[csr_matrix, csr_matrix, csr_matrix]:
    """Build sparse coefficient-to-cell-current maps (Bx, By, Bz).

    Edge selection alone determines column order; no gauge is chosen here.
    Conservation at boundaries and independence depend on that selection.
    Gather opposing-edge endpoints only for selected local edges, then assemble
    one component at a time, avoiding full (ncells, 6, 3) arrays and leaving
    mesh.cell_edge_vectors unpopulated.
    """
    _check_connectivity(mesh, topology)
    retained_edges = np.asarray(retained_edges)
    if retained_edges.ndim != 1 or (
        retained_edges.size and not np.issubdtype(retained_edges.dtype, np.integer)
    ):
        raise ValueError("retained_edges must be a one-dimensional integer array")
    retained_edges = retained_edges.astype(np.int64, copy=False)
    if np.any(retained_edges < 0) or np.any(retained_edges >= len(topology.edges)):
        raise ValueError("retained edge index is out of range")
    if len(np.unique(retained_edges)) != len(retained_edges):
        raise ValueError("retained_edges must be unique")
    ncells = len(mesh.tetrahedra)
    ndofs = len(retained_edges)
    edge_to_dof = np.full(len(topology.edges), -1, dtype=np.int64)
    edge_to_dof[retained_edges] = np.arange(ndofs)
    columns = edge_to_dof[topology.cell_edges]
    cells, local_edges = np.nonzero(columns >= 0)
    columns = columns[cells, local_edges]
    opposing_edges = np.asarray(TETRA_EDGES)[(local_edges + 3) % 6]
    tails = mesh.tetrahedra[cells, opposing_edges[:, 0]]
    heads = mesh.tetrahedra[cells, opposing_edges[:, 1]]
    scales = topology.cell_edge_signs[cells, local_edges] / (3 * mesh.volumes[cells])
    operators = []
    for component in range(3):
        values = (
            mesh.vertices[heads, component] - mesh.vertices[tails, component]
        ) * scales
        operator = coo_matrix(
            (values, (cells, columns)),
            shape=(ncells, ndofs),
        ).tocsr()
        operator.eliminate_zeros()
        operators.append(operator)
    return operators[0], operators[1], operators[2]


def reconstruct_currents(
    basis: CellCurrentBasis, coefficients: ArrayLike
) -> NDArray[np.floating | np.complexfloating]:
    """Reconstruct constant Cartesian currents from one or several coefficient sets.

    Input shape is (n_dofs,) or (n_dofs, nfields); output shape is (ncells, 3)
    or (ncells, 3, nfields), respectively. With coefficients in A and coordinates
    in metres, the result is in A/m². Complex phasors are supported.
    """
    coefficients = np.asarray(coefficients)
    if coefficients.ndim not in (1, 2) or coefficients.shape[0] != basis.n_dofs:
        raise ValueError("coefficients must have shape (n_dofs,) or (n_dofs, nfields)")
    return np.stack([operator @ coefficients for operator in basis.operators], axis=1)
