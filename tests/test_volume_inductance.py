"""Independent integral references, magnetic energy, and generalized modes."""

from math import factorial

import numpy as np
import pytest
from scipy.constants import mu_0
from scipy.linalg import eigh
from scipy.special import roots_jacobi
from scipy.sparse import csr_matrix

from bfieldtools.volume import (
    CellCurrentBasis,
    TetraMesh,
    TetraTopology,
    TreeCotreeGauge,
    build_edge_current_basis,
    inductance_matrix,
    resistance_matrix,
    tetrahedron_interaction_block,
    tetrahedron_interaction_matrix,
    tetrahedron_potential,
    tetrahedron_quadrature,
)
from bfieldtools.volume.mesh import TETRA_FACES
from .test_volume_topology import cube_mesh


def simple_tetrahedron():
    return TetraMesh(np.vstack((np.zeros(3), np.eye(3))), [[0, 1, 2, 3]])


def basis_for(mesh):
    topology = TetraTopology.from_mesh(mesh)
    return build_edge_current_basis(
        mesh, topology, TreeCotreeGauge.from_topology(topology)
    )


def triangle_rule(order):
    a, wa = roots_jacobi(order, 0, 1)
    b, wb = roots_jacobi(order, 0, 0)
    a, b = np.meshgrid((a + 1) / 2, (b + 1) / 2, indexing="ij")
    points = np.stack((1 - a, a * (1 - b), a * b), axis=-1).reshape(-1, 3)
    return points, np.outer(wa / 2, wb / 2).ravel()


def point_potential_reference(mesh, points):
    # Numerical triangle integration after integrating the radial coordinate
    # analytically. Unlike the implementation, no triangle-potential formula
    # or edge logarithm is used. Boundary faces with zero height contribute zero.
    barycentric, weights = triangle_rule(48)
    faces = mesh.vertices[mesh.tetrahedra[0, np.asarray(TETRA_FACES)]]
    area_vectors = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0]) / 2
    samples = np.einsum("qv,fvi->fqi", barycentric, faces)
    result = []
    for point in points:
        value = 0.0
        for face, area, samples_f in zip(faces, area_vectors, samples):
            factor = np.dot(face[0] - point, area)
            if abs(factor) > 1e-15:
                value += (
                    0.5
                    * factor
                    * np.sum(weights / np.linalg.norm(samples_f - point, axis=1))
                )
        result.append(value)
    return np.array(result)


def cube_self_reference():
    # Covariogram of a unit cube: integrate over displacement, then partition
    # the positive octant by its largest coordinate. The radial integral is
    # exact, leaving a smooth 2D integral independent of the tetrahedral method.
    x, w = np.polynomial.legendre.leggauss(48)
    a, b = np.meshgrid((x + 1) / 2, (x + 1) / 2, indexing="ij")
    integrand = (1 / 6 - (a + b) / 12 + a * b / 20) / np.sqrt(1 + a * a + b * b)
    return 24 * np.einsum("i,j,ij", w / 2, w / 2, integrand)


@pytest.mark.parametrize("order", [1, 2, 4])
def test_quadrature_moments(order):
    points, weights = tetrahedron_quadrature(order)
    np.testing.assert_allclose(points.sum(axis=1), 1)
    assert np.all(points > 0)
    assert np.all(weights > 0)
    np.testing.assert_allclose(weights.sum(), 1)
    for exponents in [(0, 0, 0, 0), (2 * order - 1, 0, 0, 0), (0, 1, 0, 0)]:
        exact = (
            6
            * np.prod([factorial(n) for n in exponents])
            / factorial(sum(exponents) + 3)
        )
        np.testing.assert_allclose(
            weights @ np.prod(points**exponents, axis=1), exact, rtol=1e-13
        )


def test_point_potential_at_interior_faces_edges_and_vertices():
    mesh = simple_tetrahedron()
    points = np.array(
        [[0.25, 0.25, 0.25], [0.2, 0.2, 0], [0.3, 0, 0], [0, 0, 0], [3, 2, 1]]
    )
    computed = tetrahedron_potential(mesh, points)[:, 0]
    np.testing.assert_allclose(
        computed, point_potential_reference(mesh, points), rtol=2e-10
    )
    assert np.all(np.isfinite(computed))
    assert np.all(computed > 0)
    np.testing.assert_allclose(
        tetrahedron_potential(mesh, points, point_chunk_size=1)[:, 0],
        computed,
        rtol=1e-12,
    )


def test_far_point_potential_uses_stable_regular_source_integration():
    mesh = simple_tetrahedron()
    points = np.array([[1e5, 2e5, 3e5], [5e5, 0, 0]])
    distance = np.linalg.norm(points - mesh.cell_centers[0], axis=1)
    np.testing.assert_allclose(
        tetrahedron_potential(mesh, points)[:, 0],
        mesh.volumes[0] / distance,
        rtol=1e-11,
    )


@pytest.mark.parametrize(
    "point_chunk_size,source_chunk_size", [(1, 1), (3, 2), (256, 16)]
)
def test_point_potential_batches_mixed_geometry(point_chunk_size, source_chunk_size):
    # Different scales and face normals expose incorrect pair/face packing.
    tetra = simple_tetrahedron().vertices
    vertices = np.concatenate(
        (tetra, tetra * [2, 3, 1] + [5, 0, 0], tetra[:, [1, 0, 2]] * 0.1 + 30)
    )
    mesh = TetraMesh(vertices, np.arange(12).reshape(3, 4))
    points = np.vstack((mesh.cell_centers, vertices, [[100, 200, 300]]))
    sources = np.array([2, 0, 1, 2])
    # Independent face quadrature also covers exact vertices and mixed near/far
    # classifications, without sharing the analytic triangle implementation.
    expected = np.column_stack(
        [
            point_potential_reference(
                TetraMesh(mesh.vertices[mesh.tetrahedra[cell]], [[0, 1, 2, 3]]),
                points,
            )
            for cell in sources
        ]
    )
    np.testing.assert_allclose(
        tetrahedron_potential(
            mesh,
            points,
            sources,
            point_chunk_size=point_chunk_size,
            source_chunk_size=source_chunk_size,
        ),
        expected,
        rtol=2e-8,
    )


def test_point_potential_empty_batches_and_invalid_source_chunk_size():
    mesh = simple_tetrahedron()
    assert tetrahedron_potential(mesh, np.empty((0, 3))).shape == (0, 1)
    assert tetrahedron_potential(mesh, [[0, 0, 0]], []).shape == (1, 0)
    for value in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match="source_chunk_size"):
            tetrahedron_potential(mesh, [[0, 0, 0]], source_chunk_size=value)


def test_self_integral_against_independent_face_pair_distance_reduction():
    mesh = simple_tetrahedron()
    barycentric, weights = triangle_rule(24)
    faces = mesh.vertices[mesh.tetrahedra[0, np.asarray(TETRA_FACES)]]
    area_vectors = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0]) / 2
    points = np.einsum("qv,fvi->fqi", barycentric, faces)
    reference = 0.0
    for i in range(4):
        for j in range(4):
            distances = np.linalg.norm(points[i, :, None] - points[j, None, :], axis=-1)
            reference -= (
                0.5
                * np.dot(area_vectors[i], area_vectors[j])
                * (weights @ distances @ weights)
            )
    computed = tetrahedron_interaction_matrix(mesh, rtol=1e-7, max_order=20)[0, 0]
    # The independent distance-kernel rule converges slowly on coincident faces.
    np.testing.assert_allclose(computed, reference, rtol=4e-4)


def test_cube_self_integral_against_independent_displacement_reference():
    matrix = tetrahedron_interaction_matrix(cube_mesh(1))
    np.testing.assert_allclose(matrix.sum(), cube_self_reference(), rtol=3e-6)
    np.testing.assert_array_equal(matrix, matrix.T)
    assert np.linalg.eigvalsh(matrix).min() > 0


@pytest.mark.parametrize("kind", ["face", "edge", "vertex", "near", "far"])
def test_touching_and_separated_cell_pairs(kind):
    first = simple_tetrahedron().vertices
    second = {
        "face": first * [1, 1, -1],
        "edge": first * [1, -1, -1],
        "vertex": -first,
        "near": first + [1.001, 0, 0],
        "far": first + [20, 5, 2],
    }[kind]
    mesh = TetraMesh(np.concatenate((first, second)), [[0, 1, 2, 3], [4, 5, 6, 7]])
    coarse = tetrahedron_interaction_block(mesh, [0], [1])[0, 0]
    fine = tetrahedron_interaction_block(mesh, [0], [1], rtol=1e-6, max_order=20)[0, 0]
    reverse = tetrahedron_interaction_block(mesh, [1], [0], rtol=1e-6, max_order=20)[
        0, 0
    ]
    assert fine > 0
    np.testing.assert_allclose(coarse, fine, rtol=2e-4)
    np.testing.assert_allclose(fine, reverse, rtol=1e-14)
    if kind == "far":
        barycentric, weights = tetrahedron_quadrature(5)
        cells = mesh.vertices[mesh.tetrahedra]
        p = np.einsum("qv,tvi->tqi", barycentric, cells)
        distances = np.linalg.norm(p[0, :, None] - p[1, None, :], axis=-1)
        reference = np.prod(mesh.volumes) * (weights @ (1 / distances) @ weights)
        np.testing.assert_allclose(fine, reference, rtol=1e-10)


def test_kernel_rigid_motion_and_fifth_power_scaling():
    mesh = cube_mesh(1)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    moved = TetraMesh(2 * mesh.vertices @ rotation.T + [3, 4, 5], mesh.tetrahedra)
    np.testing.assert_allclose(
        tetrahedron_interaction_matrix(moved),
        2**5 * tetrahedron_interaction_matrix(mesh),
        rtol=1e-10,
    )


def test_blocked_inductance_and_energy_against_dense_contraction():
    mesh = cube_mesh(2)
    basis = basis_for(mesh)
    kernel = tetrahedron_interaction_matrix(mesh, block_size=12)
    reference = (
        sum(operator.T @ kernel @ operator for operator in basis.operators)
        * mu_0
        / (4 * np.pi)
    )
    cached = inductance_matrix(basis, interactions=kernel, block_size=7)
    streamed = inductance_matrix(basis, block_size=12)
    np.testing.assert_allclose(cached, reference, rtol=1e-12, atol=1e-20)
    np.testing.assert_allclose(streamed, reference, rtol=1e-12, atol=1e-20)
    x = np.random.default_rng(72).normal(size=basis.n_dofs)
    currents = np.stack([operator @ x for operator in basis.operators], axis=1)
    energy = mu_0 / (8 * np.pi) * np.sum(kernel * (currents @ currents.T))
    np.testing.assert_allclose(x @ cached @ x / 2, energy)
    assert np.linalg.eigvalsh(cached).min() > 0


def test_generalized_modes_basis_invariance_and_length_scaling():
    mesh = cube_mesh(2)
    basis = basis_for(mesh)
    kernel = tetrahedron_interaction_matrix(mesh)
    magnetic = inductance_matrix(basis, interactions=kernel)
    resistance = resistance_matrix(basis, 1).toarray()
    rates, vectors = eigh(resistance, magnetic)
    assert rates.min() > 0
    np.testing.assert_allclose(
        vectors.T @ magnetic @ vectors, np.eye(basis.n_dofs), atol=2e-13
    )
    residual = resistance @ vectors - (magnetic @ vectors) * rates
    assert np.linalg.norm(residual) / np.linalg.norm(resistance @ vectors) < 1e-12
    transform, _ = np.linalg.qr(
        np.random.default_rng(6).normal(size=(basis.n_dofs, basis.n_dofs))
    )
    transformed = CellCurrentBasis(
        mesh, tuple(csr_matrix(operator @ transform) for operator in basis.operators)
    )
    np.testing.assert_allclose(
        eigh(
            resistance_matrix(transformed, 1).toarray(),
            inductance_matrix(transformed, interactions=kernel),
            eigvals_only=True,
        ),
        rates,
        rtol=1e-12,
    )
    larger = basis_for(TetraMesh(2 * mesh.vertices, mesh.tetrahedra))
    larger_magnetic = inductance_matrix(larger, interactions=32 * kernel)
    np.testing.assert_allclose(larger_magnetic, 2 * magnetic, atol=1e-20)
    np.testing.assert_allclose(
        eigh(
            resistance_matrix(larger, 1).toarray(), larger_magnetic, eigvals_only=True
        ),
        rates / 4,
        rtol=1e-12,
    )


def test_nonconvergence_is_reported_instead_of_adjusting_diagonal():
    with pytest.raises(RuntimeError, match="did not converge"):
        tetrahedron_interaction_matrix(
            simple_tetrahedron(), quadrature_order=1, max_order=3, rtol=1e-12
        )


@pytest.mark.parametrize(
    "options",
    [
        {"quadrature_order": 0},
        {"quadrature_order": 1.5},
        {"rtol": 0},
        {"rtol": np.nan},
        {"max_order": 4},
        {"block_size": 0},
    ],
)
def test_invalid_integration_controls(options):
    with pytest.raises(ValueError):
        tetrahedron_interaction_matrix(simple_tetrahedron(), **options)


@pytest.mark.parametrize(
    "kernel", [np.ones((2, 2)), np.array([[np.nan]]), np.array([[1j]])]
)
def test_invalid_precomputed_kernel(kernel):
    with pytest.raises(ValueError, match="interactions"):
        inductance_matrix(basis_for(simple_tetrahedron()), interactions=kernel)


def test_empty_current_space():
    assert inductance_matrix(basis_for(simple_tetrahedron())).shape == (0, 0)
