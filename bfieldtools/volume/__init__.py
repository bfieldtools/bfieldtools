"""Tetrahedral meshes and tools for volumetric conductor models."""

from .mesh import TetraMesh, tetrahedralize, load_example_tetra_mesh
from .dynamics import loop_coupling, modal_response
from .currents import ConservedFluxBasis, augmented_flux_basis, dual_cycle_basis
from .basis import (
    CellCurrentBasis,
    build_edge_current_basis,
    build_face_current_basis,
    cell_current_operators,
    cell_edge_current_vectors,
    reconstruct_currents,
)
from .impedance import inductance_matrix, resistance_matrix
from .integrals import (
    iter_tetrahedron_interactions,
    tetrahedron_interaction_block,
    tetrahedron_interaction_matrix,
    tetrahedron_potential,
    tetrahedron_quadrature,
)
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
    "load_example_tetra_mesh",
    "loop_coupling",
    "modal_response",
    "TetraTopology",
    "TreeCotreeGauge",
    "build_topology",
    "gradient_matrix",
    "curl_matrix",
    "divergence_matrix",
    "tree_cotree_gauge",
    "CellCurrentBasis",
    "ConservedFluxBasis",
    "dual_cycle_basis",
    "augmented_flux_basis",
    "build_face_current_basis",
    "build_edge_current_basis",
    "cell_current_operators",
    "cell_edge_current_vectors",
    "reconstruct_currents",
    "resistance_matrix",
    "inductance_matrix",
    "tetrahedron_potential",
    "tetrahedron_quadrature",
    "tetrahedron_interaction_block",
    "tetrahedron_interaction_matrix",
    "iter_tetrahedron_interactions",
]
