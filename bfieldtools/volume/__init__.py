"""Tetrahedral meshes and tools for volumetric conductor models."""

from .mesh import TetraMesh, tetrahedralize
from .basis import (
    CellCurrentBasis,
    build_edge_current_basis,
    cell_current_operators,
    cell_edge_current_vectors,
    reconstruct_currents,
)
from .impedance import resistance_matrix
from .topology import (
    TetraTopology,
    TreeCotreeGauge,
    build_topology,
    curl_matrix,
    divergence_matrix,
    gradient_matrix,
    tree_cotree_gauge,
)

__all__ = [
    "TetraMesh",
    "tetrahedralize",
    "TetraTopology",
    "TreeCotreeGauge",
    "build_topology",
    "gradient_matrix",
    "curl_matrix",
    "divergence_matrix",
    "tree_cotree_gauge",
    "CellCurrentBasis",
    "build_edge_current_basis",
    "cell_current_operators",
    "cell_edge_current_vectors",
    "reconstruct_currents",
    "resistance_matrix",
]
