"""Physical orientation, conservation, and independence of volume edge loops."""

from itertools import permutations, product

import numpy as np
import pytest

from bfieldtools.volume import TetraMesh, build_topology, tree_cotree_gauge
from bfieldtools.volume.mesh import TETRA_EDGES


def cube_mesh(n=2, omit=()):
    """Conforming Freudenthal subdivision; optionally remove whole voxels."""
    shape = (n + 1,) * 3
    vertices = np.array(list(product(range(n + 1), repeat=3)), dtype=float)
    cells = []
    for origin in product(range(n), repeat=3):
        if origin in omit:
            continue
        for axes in permutations(range(3)):
            points = [np.array(origin)]
            for axis in axes:
                point = points[-1].copy()
                point[axis] += 1
                points.append(point)
            cells.append([np.ravel_multi_index(point, shape) for point in points])
    return TetraMesh(vertices / n, cells)


def assert_zero(matrix):
    matrix = matrix.copy()
    matrix.eliminate_zeros()
    assert matrix.nnz == 0


def test_incidence_and_geometry():
    mesh = cube_mesh()
    topology = build_topology(mesh)
    assert_zero(topology.curl @ topology.gradient)
    assert_zero(topology.divergence @ topology.curl)

    np.testing.assert_allclose(
        topology.gradient @ mesh.vertices,
        mesh.vertices[topology.edges[:, 1]] - mesh.vertices[topology.edges[:, 0]],
    )
    local = mesh.tetrahedra[:, TETRA_EDGES]
    np.testing.assert_array_equal(
        np.sort(local, axis=2), topology.edges[topology.cell_edges]
    )
    d = mesh.cell_edge_vectors
    global_d = topology.gradient @ mesh.vertices
    np.testing.assert_allclose(
        d, topology.cell_edge_signs[:, :, None] * global_d[topology.cell_edges]
    )
    np.testing.assert_array_equal(
        topology.faces[topology.cell_faces],
        np.sort(
            mesh.tetrahedra[:, ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))], axis=2
        ),
    )

    # Independent geometric flux calculation, using globally oriented normals.
    # This checks the opposing-edge curls and BOTH incidence orientation maps.
    for cell, nodes in enumerate(mesh.tetrahedra):
        points = mesh.vertices[nodes]
        gradients = np.empty((4, 3))
        gradients[1:] = np.linalg.inv((points[1:] - points[0]).T)
        gradients[0] = -gradients[1:].sum(axis=0)
        edges = np.asarray(TETRA_EDGES)
        curls = 2 * np.cross(gradients[edges[:, 0]], gradients[edges[:, 1]])
        np.testing.assert_allclose(
            curls, d[cell, (np.arange(6) + 3) % 6] / (3 * mesh.volumes[cell])
        )
        face_ids = topology.cell_faces[cell]
        triangles = mesh.vertices[topology.faces[face_ids]]
        area_vectors = (
            np.cross(
                triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
            )
            / 2
        )
        measured = area_vectors @ (curls * topology.cell_edge_signs[cell, :, None]).T
        np.testing.assert_allclose(
            measured,
            topology.curl[face_ids][:, topology.cell_edges[cell]].toarray(),
            atol=1e-14,
        )


@pytest.mark.parametrize("n", [1, 2, 3])
def test_gauge_is_complete_and_independent(n):
    mesh = cube_mesh(n)
    topology = build_topology(mesh)
    gauge = tree_cotree_gauge(topology)
    interior_vertices = topology.active_vertex_mask & ~topology.boundary_vertex_mask
    assert len(gauge.tree_edges) == np.count_nonzero(interior_vertices)
    assert len(gauge.retained_edges) == np.count_nonzero(
        ~topology.boundary_edge_mask
    ) - np.count_nonzero(interior_vertices)
    assert not np.any(topology.boundary_edge_mask[gauge.tree_edges])
    assert_zero(gauge.curl[topology.boundary_face_mask])
    assert_zero(topology.divergence @ gauge.curl)
    assert_zero(topology.curl @ gauge.prolongation - gauge.curl)
    assert np.linalg.matrix_rank(gauge.curl.toarray()) == gauge.n_dofs

    # Compare with ALL conserved face-flux fields, independently of edge curls.
    Di = topology.divergence[:, ~topology.boundary_face_mask].toarray()
    assert gauge.n_dofs == Di.shape[1] - np.linalg.matrix_rank(Di)
    random = np.random.default_rng(17).normal(size=Di.shape[1])
    conserved = random - np.linalg.lstsq(Di, Di @ random, rcond=None)[0]
    reconstructed = np.linalg.lstsq(
        gauge.curl[~topology.boundary_face_mask].toarray(), conserved, rcond=None
    )[0]
    np.testing.assert_allclose(
        gauge.curl[~topology.boundary_face_mask] @ reconstructed, conserved, atol=1e-12
    )


def test_boundary_endpoints_do_not_make_an_edge_a_boundary_edge():
    topology = build_topology(cube_mesh(1))
    gauge = tree_cotree_gauge(topology)
    # The body diagonal is interior even though both endpoints are on the surface.
    assert gauge.n_dofs == 1
    assert len(gauge.tree_edges) == 0
    assert np.all(topology.boundary_vertex_mask[topology.edges[gauge.retained_edges]])


def test_disconnected_components_and_unused_vertex():
    first = cube_mesh()
    mesh = TetraMesh(
        np.concatenate((first.vertices, first.vertices + [2, 0, 0], [[9, 9, 9]])),
        np.concatenate((first.tetrahedra, first.tetrahedra + len(first.vertices))),
    )
    topology = build_topology(mesh)
    gauge = tree_cotree_gauge(topology)
    assert topology.n_components == 2
    assert topology.vertex_components[-1] == -1
    assert len(gauge.tree_edges) == 2
    assert gauge.n_dofs == 2 * tree_cotree_gauge(build_topology(first)).n_dofs
    assert_zero(topology.divergence @ gauge.curl)
    assert_zero(gauge.curl[topology.boundary_face_mask])


def test_single_tetrahedron_has_no_isolated_current():
    mesh = TetraMesh([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 2, 3]])
    topology = build_topology(mesh)
    gauge = tree_cotree_gauge(topology)
    assert gauge.n_dofs == 0
    assert gauge.prolongation.shape == (6, 0)
    assert gauge.curl.shape == (4, 0)


def test_handle_needs_global_patterns():
    # A square ring extruded through the cube: one connected genus-one boundary.
    omit = tuple((1, 1, z) for z in range(3))
    topology = build_topology(cube_mesh(3, omit=omit))
    assert_zero(topology.divergence @ topology.curl)
    with pytest.raises(ValueError, match="handles"):
        tree_cotree_gauge(topology)


def test_cavity_needs_additional_gauge_constraints():
    topology = build_topology(cube_mesh(3, omit=((1, 1, 1),)))
    with pytest.raises(ValueError, match="cavities"):
        tree_cotree_gauge(topology)


def test_nonmanifold_boundary_is_rejected():
    # Two tetrahedra touch only along an edge; four boundary faces meet there.
    mesh = TetraMesh(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, -1, 0], [0, 0, -1]],
        [[0, 1, 2, 3], [0, 1, 4, 5]],
    )
    with pytest.raises(ValueError, match="manifold"):
        tree_cotree_gauge(build_topology(mesh))


def test_geometry_properties_are_cached_and_readonly():
    mesh = cube_mesh()
    for name in (
        "cell_centers",
        "cell_radii",
        "cell_edge_vectors",
        "cell_quality",
        "spatial_cell_order",
    ):
        array = getattr(mesh, name)
        assert array is getattr(mesh, name)
        assert not array.flags.writeable


def test_incidence_builders_are_separate_and_cached():
    from bfieldtools.volume import gradient_matrix, curl_matrix, divergence_matrix

    topology = build_topology(cube_mesh())
    for name, builder in (
        ("gradient", gradient_matrix),
        ("curl", curl_matrix),
        ("divergence", divergence_matrix),
    ):
        assert name not in topology.__dict__
        explicit = builder(topology)
        cached = getattr(topology, name)
        assert_zero(explicit - cached)
        assert cached is getattr(topology, name)


def test_dataclass_construction_helpers():
    from dataclasses import is_dataclass
    from bfieldtools.volume import TetraTopology, TreeCotreeGauge

    mesh = cube_mesh()
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    reference = tree_cotree_gauge(build_topology(mesh))
    assert is_dataclass(topology)
    assert is_dataclass(gauge)
    np.testing.assert_array_equal(gauge.retained_edges, reference.retained_edges)
    np.testing.assert_array_equal(gauge.tree_edges, reference.tree_edges)
    assert_zero(gauge.curl - reference.curl)
