"""Complete volume currents on handles, cavities, and disconnected solids."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import eigvalsh, null_space
from scipy.sparse import csr_matrix

from bfieldtools.volume import (
    ConservedFluxBasis,
    TetraMesh,
    augmented_flux_basis,
    build_edge_current_basis,
    build_face_current_basis,
    build_topology,
    inductance_matrix,
    loop_coupling,
    reconstruct_currents,
    resistance_matrix,
    tetrahedron_interaction_matrix,
    tree_cotree_gauge,
)
from .test_volume_topology import cube_mesh


def ring_mesh():
    return cube_mesh(3, omit=tuple((1, 1, z) for z in range(3)))


def double_mesh(mesh):
    return TetraMesh(
        np.concatenate((mesh.vertices, mesh.vertices + [2, 0, 0], [[9, 9, 9]])),
        np.concatenate((mesh.tetrahedra, mesh.tetrahedra + len(mesh.vertices))),
    )


def conserved_reference(topology):
    """Independent full divergence-nullspace reference; no local curl input."""
    interior = np.flatnonzero(~topology.boundary_face_mask)
    z = null_space(topology.divergence[:, interior].toarray().astype(float))
    fluxes = np.zeros((len(topology.faces), z.shape[1]))
    fluxes[interior] = z
    return ConservedFluxBasis(topology, csr_matrix(fluxes))


@pytest.mark.parametrize(
    "mesh, expected_global",
    [
        (cube_mesh(1), 0),
        (cube_mesh(2), 0),
        (ring_mesh(), 1),
        (cube_mesh(3, omit=((1, 1, 1),)), 0),
        (double_mesh(ring_mesh()), 2),
        (cube_mesh(5, omit=tuple((x, 2, z) for x in (1, 3) for z in range(5))), 2),
        (cube_mesh(5, omit=tuple((3, 3, z) for z in range(5)) + ((1, 1, 2),)), 1),
        # Thin two-handle slab: no interior vertices, so the stack is wide.
        (
            cube_mesh(
                5,
                omit=tuple(
                    (x, y, z)
                    for x in range(5)
                    for y in range(5)
                    for z in range(5)
                    if z != 0 or (x, y) in ((1, 2), (3, 2))
                ),
            ),
            2,
        ),
        (TetraMesh(np.vstack((np.zeros(3), np.eye(3))), [[0, 1, 2, 3]]), 0),
    ],
)
def test_complete_flux_spaces(mesh, expected_global):
    topology = build_topology(mesh)
    augmented = augmented_flux_basis(topology)
    assert augmented.n_global == expected_global
    reference = conserved_reference(topology)
    dimension = reference.n_dofs
    assert augmented.n_dofs == dimension
    fluxes = augmented.fluxes.toarray()
    if dimension:
        assert np.linalg.matrix_rank(fluxes) == dimension
    else:
        assert fluxes.shape[1] == 0
    np.testing.assert_allclose(topology.divergence @ fluxes, 0, atol=2e-13)
    np.testing.assert_array_equal(fluxes[topology.boundary_face_mask], 0)
    coefficients = np.linalg.lstsq(fluxes, reference.fluxes.toarray(), rcond=None)[0]
    np.testing.assert_allclose(
        fluxes @ coefficients, reference.fluxes.toarray(), atol=2e-13
    )
    # Independently integrate reconstructed currents through every face of
    # every cell: checks RT reconstruction, orientation, and continuity.
    basis = build_face_current_basis(mesh, augmented)
    currents = reconstruct_currents(basis, np.eye(dimension))
    triangles = mesh.vertices[topology.faces]
    area_vectors = (
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        / 2
    )
    measured = np.einsum("tfi,tik->tfk", area_vectors[topology.cell_faces], currents)
    np.testing.assert_allclose(measured, fluxes[topology.cell_faces], atol=3e-13)
    nlocal = len(augmented.local_edges)
    np.testing.assert_array_equal(
        fluxes[:, :nlocal], topology.curl[:, augmented.local_edges].toarray()
    )
    global_fluxes = fluxes[:, nlocal:]
    np.testing.assert_allclose(
        global_fluxes.T @ global_fluxes, np.eye(expected_global), atol=1e-13
    )
    # The stack uses ALL interior edge curls, including redundant potentials.
    np.testing.assert_allclose(
        topology.curl[:, ~topology.boundary_edge_mask].T @ global_fluxes,
        0,
        atol=1e-13,
    )


def test_ring_global_current_crosses_cut_local_curls_do_not():
    mesh = ring_mesh()
    topology = build_topology(mesh)
    space = augmented_flux_basis(topology)
    triangles = mesh.vertices[topology.faces]
    on_cut = np.all(np.isclose(triangles[:, :, 1], 1 / 3), axis=1)
    on_cut &= np.all(triangles[:, :, 0] >= 2 / 3 - 1e-14, axis=1)
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    cut = on_cut * np.sign(normals[:, 1])
    periods = cut @ space.fluxes
    np.testing.assert_allclose(periods[:-1], 0, atol=1e-14)
    assert abs(periods[-1]) > 0.1
    currents = reconstruct_currents(
        build_face_current_basis(mesh, space),
        np.r_[np.zeros(space.n_dofs - 1), 1 / periods[-1]],
    )
    assert np.linalg.norm(currents) > 0
    assert (cut @ (space.fluxes[:, -1] / periods[-1]))[0] == pytest.approx(1)


def test_face_reconstruction_matches_legacy_edge_operators():
    mesh = cube_mesh(2)
    topology = build_topology(mesh)
    legacy = build_edge_current_basis(mesh, topology, tree_cotree_gauge(topology))
    # Use the old gauge's exact flux columns, testing the new reconstruction
    # without a numerically fitted change of basis.
    flux_basis = ConservedFluxBasis(topology, tree_cotree_gauge(topology).curl)
    reconstructed = build_face_current_basis(mesh, flux_basis)
    for actual, expected in zip(reconstructed.operators, legacy.operators):
        np.testing.assert_allclose(actual.toarray(), expected.toarray(), atol=2e-14)


def test_impedance_coupling_and_modes_are_basis_invariant():
    mesh = ring_mesh()
    topology = build_topology(mesh)
    reference, augmented = conserved_reference(topology), augmented_flux_basis(topology)
    independent = build_face_current_basis(mesh, reference)
    local_global = build_face_current_basis(mesh, augmented)
    transform = np.linalg.lstsq(
        reference.fluxes.toarray(), augmented.fluxes.toarray(), rcond=None
    )[0]
    rho = np.linspace(1, 2, len(mesh.tetrahedra))
    rd = resistance_matrix(independent, rho).toarray()
    ra = resistance_matrix(local_global, rho).toarray()
    np.testing.assert_allclose(ra, transform.T @ rd @ transform, atol=3e-12)
    kernel = tetrahedron_interaction_matrix(mesh, quadrature_order=2, touching_order=4)
    md = inductance_matrix(independent, interactions=kernel)
    ma = inductance_matrix(local_global, interactions=kernel)
    np.testing.assert_allclose(ma, transform.T @ md @ transform, atol=1e-18)
    assert np.linalg.eigvalsh(ra).min() > 0
    np.testing.assert_allclose(eigvalsh(ra, ma), eigvalsh(rd, md), rtol=1e-10)
    loop = [[-1, -1, 2], [2, -1, 2], [2, 2, 2], [-1, 2, 2], [-1, -1, 2]]
    np.testing.assert_allclose(
        loop_coupling(local_global, loop),
        transform.T @ loop_coupling(independent, loop),
        atol=1e-19,
    )
    # Generic coefficients include the global mode and all local/global terms.
    x = np.random.default_rng(2).normal(size=local_global.n_dofs) * (1 + 1j)
    currents = reconstruct_currents(local_global, x)
    np.testing.assert_allclose(
        np.vdot(x, ra @ x).real,
        np.sum(rho[:, None] * mesh.volumes[:, None] * abs(currents) ** 2),
    )


def test_orientation_permutation_and_geometric_scaling():
    mesh = ring_mesh()
    space = augmented_flux_basis(build_topology(mesh))
    original = build_face_current_basis(mesh, space)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    moved = TetraMesh(2 * mesh.vertices @ rotation.T + [3, 4, 5], mesh.tetrahedra)
    scaled = build_face_current_basis(moved, space)
    x = np.random.default_rng(3).normal(size=space.n_dofs)
    np.testing.assert_allclose(
        reconstruct_currents(scaled, x),
        reconstruct_currents(original, x) @ rotation.T / 4,
        atol=1e-13,
    )
    # Reversing cell vertex order is repaired by TetraMesh; arbitrary vertex
    # and cell permutations must preserve the physical represented space.
    rng = np.random.default_rng(9)
    order = rng.permutation(len(mesh.vertices))
    inverse = np.argsort(order)
    cell_order = rng.permutation(len(mesh.tetrahedra))
    shuffled = TetraMesh(
        mesh.vertices[order], inverse[mesh.tetrahedra[cell_order, ::-1]]
    )
    shuffled_space = augmented_flux_basis(build_topology(shuffled))
    assert shuffled_space.n_global == 1
    shuffled_basis = build_face_current_basis(shuffled, shuffled_space)
    a = np.vstack([operator.toarray() for operator in original.operators])
    b = np.vstack(
        [
            operator.toarray()[np.argsort(cell_order)]
            for operator in shuffled_basis.operators
        ]
    )
    fit = np.linalg.lstsq(b, a, rcond=None)[0]
    np.testing.assert_allclose(b @ fit, a, atol=2e-12)


def test_reject_nonmanifold_boundaries():
    # Edge pinch, then a vertex-only pinch with otherwise closed edge links.
    vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, -1, 0], [0, 0, -1]]
    mesh = TetraMesh(vertices, [[0, 1, 2, 3], [0, 1, 4, 5]])
    with pytest.raises(ValueError, match="manifold"):
        augmented_flux_basis(build_topology(mesh))
    vertices.append([-1, 0, 0])
    mesh = TetraMesh(vertices, [[0, 1, 2, 3], [0, 6, 4, 5]])
    with pytest.raises(ValueError, match="manifold"):
        augmented_flux_basis(build_topology(mesh))


def test_dense_guard_and_invalid_tolerance():
    topology = build_topology(ring_mesh())
    with pytest.raises(ValueError, match="max_dense_entries"):
        augmented_flux_basis(topology, max_dense_entries=1)
    for rtol in (0, -1, 1, np.nan):
        with pytest.raises(ValueError, match="rtol"):
            augmented_flux_basis(topology, rtol=rtol)
    for size in (0, 1.5, True):
        with pytest.raises(ValueError, match="max_dense_entries"):
            augmented_flux_basis(topology, max_dense_entries=size)


def test_reconstruction_rejects_invalid_fluxes_and_connectivity():
    mesh = cube_mesh(1)
    topology = build_topology(mesh)
    space = augmented_flux_basis(topology)
    bad = space.fluxes.astype(float).tolil()
    bad[np.flatnonzero(topology.boundary_face_mask)[0], 0] = 1
    with pytest.raises(ValueError, match="boundary"):
        build_face_current_basis(mesh, replace(space, fluxes=bad.tocsr()))
    bad = space.fluxes.astype(float).tolil()
    bad[np.flatnonzero(~topology.boundary_face_mask)[0], 0] += 1
    with pytest.raises(ValueError, match="conserve"):
        build_face_current_basis(mesh, replace(space, fluxes=bad.tocsr()))
    bad = space.fluxes.astype(float)
    bad.data[0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        build_face_current_basis(mesh, replace(space, fluxes=bad))
    with pytest.raises(ValueError, match="one row per face"):
        build_face_current_basis(mesh, replace(space, fluxes=csr_matrix((1, 1))))
    shuffled = TetraMesh(mesh.vertices, mesh.tetrahedra[::-1])
    with pytest.raises(ValueError, match="connectivity"):
        build_face_current_basis(shuffled, space)
