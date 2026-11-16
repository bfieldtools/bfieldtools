"""Independent checks of volume-current reconstruction and Joule-loss assembly."""

import numpy as np
import pytest
from scipy.sparse import csr_matrix, isspmatrix_csr

from bfieldtools.volume import (
    CellCurrentBasis,
    TreeCotreeGauge,
    build_edge_current_basis,
    cell_current_operators,
    TetraMesh,
    TetraTopology,
    reconstruct_currents,
    resistance_matrix,
)
from bfieldtools.volume.mesh import TETRA_EDGES
from .test_volume_topology import cube_mesh


def make_basis(mesh):
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    return build_edge_current_basis(mesh, topology, gauge)


def barycentric_currents(mesh, topology, gauge, coefficients):
    """Independent gradient-based reconstruction, without opposing-edge vectors."""
    t = gauge.prolongation @ coefficients
    result = []
    for cell, nodes in enumerate(mesh.tetrahedra):
        points = mesh.vertices[nodes]
        gradients = np.empty((4, 3))
        gradients[1:] = np.linalg.inv((points[1:] - points[0]).T)
        gradients[0] = -gradients[1:].sum(axis=0)
        local = t[topology.cell_edges[cell]] * topology.cell_edge_signs[cell]
        result.append(
            sum(
                2 * local[n] * np.cross(gradients[i], gradients[j])
                for n, (i, j) in enumerate(TETRA_EDGES)
            )
        )
    return np.array(result)


def test_analytic_cube_body_diagonal_resistance():
    # Six equal tetrahedra, one retained edge (the cube body diagonal).
    # For unit coefficient, each cell has |J|=2, V=1/6, hence R=4*rho.
    basis = make_basis(cube_mesh(1))
    matrix = resistance_matrix(basis, 2.5)
    assert isspmatrix_csr(matrix)
    np.testing.assert_allclose(matrix.toarray(), [[10.0]])


def test_joule_loss_and_bilinear_form_against_independent_currents():
    basis = make_basis(cube_mesh())
    mesh = basis.mesh
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    rng = np.random.default_rng(51)
    rho = rng.uniform(0.1, 3, len(mesh.tetrahedra))
    matrix = resistance_matrix(basis, rho)
    # Assembly must leave the full six-edge geometry cache unpopulated.
    assert "cell_edge_vectors" not in mesh.__dict__
    x, y = rng.normal(size=(2, basis.n_dofs))
    jx = barycentric_currents(mesh, topology, gauge, x)
    jy = barycentric_currents(mesh, topology, gauge, y)
    np.testing.assert_allclose(reconstruct_currents(basis, x), jx, atol=1e-13)
    np.testing.assert_allclose(
        x @ matrix @ x, np.sum(rho * mesh.volumes * np.sum(jx * jx, axis=1))
    )
    np.testing.assert_allclose(
        x @ matrix @ y, np.sum(rho * mesh.volumes * np.sum(jx * jy, axis=1))
    )
    np.testing.assert_allclose(matrix.toarray(), matrix.toarray().T, atol=1e-14)
    assert np.linalg.eigvalsh(matrix.toarray()).min() > 0


def test_reconstructed_flux_matches_incidence_on_every_cell_face():
    basis = make_basis(cube_mesh())
    topology = TetraTopology.from_mesh(basis.mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    coefficients = np.random.default_rng(3).normal(size=basis.n_dofs)
    currents = reconstruct_currents(basis, coefficients)
    fluxes = gauge.curl @ coefficients
    triangles = basis.mesh.vertices[topology.faces]
    areas = (
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        / 2
    )
    measured = np.einsum("tfi,ti->tf", areas[topology.cell_faces], currents)
    np.testing.assert_allclose(measured, fluxes[topology.cell_faces], atol=1e-13)
    np.testing.assert_array_equal(fluxes[topology.boundary_face_mask], 0)
    np.testing.assert_allclose(topology.divergence @ fluxes, 0, atol=1e-13)


def test_multiple_fields_and_complex_phasors():
    basis = make_basis(cube_mesh())
    rng = np.random.default_rng(14)
    coefficients = rng.normal(size=(basis.n_dofs, 3)) + 1j * rng.normal(
        size=(basis.n_dofs, 3)
    )
    currents = reconstruct_currents(basis, coefficients)
    assert currents.shape == (len(basis.mesh.tetrahedra), 3, 3)
    for i in range(3):
        np.testing.assert_allclose(
            currents[:, :, i], reconstruct_currents(basis, coefficients[:, i])
        )
    matrix = resistance_matrix(basis, 1.0)
    x = coefficients[:, 0]
    np.testing.assert_allclose(
        np.vdot(x, matrix @ x).real,
        np.sum(basis.mesh.volumes[:, None] * abs(currents[:, :, 0]) ** 2),
    )


def test_material_scaling_linearity_and_geometric_scaling():
    mesh = cube_mesh()
    basis = make_basis(mesh)
    matrix = resistance_matrix(basis, 1).toarray()
    np.testing.assert_allclose(resistance_matrix(basis, 3).toarray(), 3 * matrix)
    np.testing.assert_allclose(
        resistance_matrix(basis, np.ones(len(mesh.tetrahedra))).toarray(), matrix
    )
    rho = np.linspace(1, 2, len(mesh.tetrahedra))
    np.testing.assert_allclose(
        resistance_matrix(basis, rho + 1).toarray(),
        resistance_matrix(basis, rho).toarray() + matrix,
    )
    # Proper rotation, translation, and uniform scaling preserve the same edge
    # coefficients/tree. R scales as inverse length, and J as inverse area.
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    transformed = make_basis(
        TetraMesh(2 * mesh.vertices @ rotation.T + [3, 4, 5], mesh.tetrahedra)
    )
    np.testing.assert_allclose(
        resistance_matrix(transformed, 1).toarray(), matrix / 2, atol=1e-14
    )
    x = np.random.default_rng(7).normal(size=basis.n_dofs)
    np.testing.assert_allclose(
        reconstruct_currents(transformed, x),
        reconstruct_currents(basis, x) @ rotation.T / 4,
        atol=1e-13,
    )


@pytest.mark.parametrize("rho", [0, -1, np.nan, np.inf, 1j, [1], [[1]], [1, -1]])
def test_invalid_resistivity(rho):
    basis = make_basis(cube_mesh())
    with pytest.raises(ValueError, match="resistivity"):
        resistance_matrix(basis, rho)


@pytest.mark.parametrize("coefficients", [0.0, [1.0], np.zeros((2, 3, 4))])
def test_invalid_coefficients(coefficients):
    with pytest.raises(ValueError, match="coefficients"):
        reconstruct_currents(make_basis(cube_mesh()), coefficients)


def test_empty_isolated_current_space():
    mesh = TetraMesh([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 2, 3]])
    basis = make_basis(mesh)
    assert resistance_matrix(basis, 1).shape == (0, 0)
    np.testing.assert_array_equal(
        reconstruct_currents(basis, np.empty(0)), np.zeros((1, 3))
    )


def test_general_cell_basis_has_the_same_physical_bilinear_form():
    basis = make_basis(cube_mesh())
    rng = np.random.default_rng(71)
    transform = rng.normal(size=(basis.n_dofs, 4))
    mixed = CellCurrentBasis(
        basis.mesh,
        tuple(csr_matrix(operator @ transform) for operator in basis.operators),
    )
    coefficients = rng.normal(size=4)
    np.testing.assert_allclose(
        reconstruct_currents(mixed, coefficients),
        reconstruct_currents(basis, transform @ coefficients),
        atol=1e-13,
    )
    np.testing.assert_allclose(
        resistance_matrix(mixed, 2).toarray(),
        transform.T @ resistance_matrix(basis, 2) @ transform,
        atol=1e-12,
    )


def test_topology_and_gauge_reused_with_new_geometry():
    mesh = cube_mesh()
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    basis = build_edge_current_basis(mesh, topology, gauge)
    scaled = TetraMesh(2 * mesh.vertices, mesh.tetrahedra)
    scaled_basis = build_edge_current_basis(scaled, topology, gauge)
    np.testing.assert_allclose(
        resistance_matrix(scaled_basis, 1).toarray(),
        resistance_matrix(basis, 1).toarray() / 2,
        atol=1e-14,
    )


def test_edge_selection_controls_column_order():
    mesh = cube_mesh()
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    basis = build_edge_current_basis(mesh, topology, gauge)
    permutation = np.random.default_rng(9).permutation(basis.n_dofs)
    reordered = cell_current_operators(
        mesh, topology, gauge.retained_edges[permutation]
    )
    for expected, actual in zip(basis.operators, reordered):
        np.testing.assert_allclose(actual.toarray(), expected[:, permutation].toarray())


def test_geometry_and_topology_cell_order_must_match():
    mesh = cube_mesh()
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    reordered = TetraMesh(mesh.vertices, mesh.tetrahedra[::-1])
    with pytest.raises(ValueError, match="ordered connectivity"):
        build_edge_current_basis(reordered, topology, gauge)


@pytest.mark.parametrize("selection", [[-1], [99999], [1, 1], [[1]], [1.5]])
def test_invalid_edge_selection(selection):
    mesh = cube_mesh()
    with pytest.raises(ValueError, match="edge"):
        cell_current_operators(mesh, TetraTopology.from_mesh(mesh), selection)


@pytest.mark.parametrize("case", ["count", "rows", "columns", "complex", "nonfinite"])
def test_cell_basis_checks_operator_layout_and_values(case):
    mesh = cube_mesh(1)
    operators = [csr_matrix(np.ones((len(mesh.tetrahedra), 1))) for _ in range(3)]
    if case == "count":
        operators.pop()
    elif case == "rows":
        operators = [operator[:-1] for operator in operators]
    elif case == "columns":
        operators[0] = csr_matrix((len(mesh.tetrahedra), 2))
    elif case == "complex":
        operators[0] = operators[0].astype(complex)
    else:
        operators[0].data[0] = np.nan
    with pytest.raises(ValueError, match="operator"):
        CellCurrentBasis(mesh, tuple(operators))
