"""Constant cell currents from gauge-fixed Whitney edge potentials."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import coo_matrix, csr_matrix

from .mesh import TETRA_EDGES
from .topology import TetraTopology, TreeCotreeGauge, tree_cotree_gauge


@dataclass(frozen=True, eq=False)
class EdgeCurrentBasis:
    """Independent edge-current basis and cached component operators.

    Construct with from_topology(). Coefficients are the retained, globally
    oriented edge integrals, in gauge.retained_edges order, with units A.
    Each operator maps these coefficients to one Cartesian current-density
    component in every cell. Geometry must use metres for SI results.
    """

    topology: TetraTopology
    gauge: TreeCotreeGauge

    @classmethod
    def from_topology(cls, topology: TetraTopology) -> EdgeCurrentBasis:
        """Build the boundary-compatible tree gauge through a separate function."""
        return build_edge_current_basis(topology)

    @cached_property
    def operators(self) -> tuple[csr_matrix, csr_matrix, csr_matrix]:
        """Tuple (Bx, By, Bz) of CSR matrices, each (ncells, n_dofs)."""
        return cell_current_operators(self)


def build_edge_current_basis(topology: TetraTopology) -> EdgeCurrentBasis:
    """Choose the independent solid-ball edge basis for a mesh topology."""
    return EdgeCurrentBasis(topology, tree_cotree_gauge(topology))


def cell_edge_current_vectors(topology: TetraTopology) -> NDArray[np.float64]:
    """Cartesian edge-current vectors in each cell, shape (ncells, 6, 3).

    Entry [T,n] is the curl of the GLOBAL function for topology.cell_edges[T,n].
    With local ordering (01,02,03,23,31,12), it is the opposing directed edge
    divided by 3V, multiplied by the potential edge's local-to-global sign.
    """
    mesh = topology.mesh
    return (
        mesh.cell_edge_vectors[:, (np.arange(6) + 3) % 6]
        * topology.cell_edge_signs[:, :, None]
        / (3 * mesh.volumes[:, None, None])
    )


def cell_current_operators(
    basis: EdgeCurrentBasis,
) -> tuple[csr_matrix, csr_matrix, csr_matrix]:
    """Build sparse coefficient-to-cell-current maps (Bx, By, Bz).

    Gather opposing-edge endpoints only for retained local edges, then assemble
    one Cartesian component at a time. This avoids full (ncells, 6, 3) vector
    arrays and does not populate mesh.cell_edge_vectors. Boundary and tree
    coefficients are identically zero.
    """
    topology = basis.topology
    mesh = topology.mesh
    ncells = len(mesh.tetrahedra)
    ndofs = basis.gauge.n_dofs
    edge_to_dof = np.full(len(topology.edges), -1, dtype=np.int64)
    edge_to_dof[basis.gauge.retained_edges] = np.arange(ndofs)
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
    basis: EdgeCurrentBasis, coefficients: ArrayLike
) -> NDArray[np.floating | np.complexfloating]:
    """Reconstruct constant Cartesian currents from one or several coefficient sets.

    Input shape is (n_dofs,) or (n_dofs, nfields); output shape is (ncells, 3)
    or (ncells, 3, nfields), respectively. With coefficients in A and coordinates
    in metres, the result is in A/m². Complex phasors are supported.
    """
    coefficients = np.asarray(coefficients)
    if coefficients.ndim not in (1, 2) or coefficients.shape[0] != basis.gauge.n_dofs:
        raise ValueError("coefficients must have shape (n_dofs,) or (n_dofs, nfields)")
    return np.stack([operator @ coefficients for operator in basis.operators], axis=1)
