"""Numerical agreement and cache lifetime for the cached integral interface."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from bfieldtools import integrals
from bfieldtools import integrals_cached
from bfieldtools.integrals_cached import TriangleIntegrals


@pytest.fixture
def geometry():
    triangles = np.array(
        [
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            [[2, 0, 0], [4, 0, 0], [2, 3, 0]],
            [[0, 2, 1], [1, 2, 2], [0, 4, 1]],
        ],
        dtype=float,
    )
    area_vectors = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    areas = np.linalg.norm(area_vectors, axis=-1) / 2
    normals = area_vectors / (2 * areas[:, None])
    points = np.vstack(([[0.2, 0.3, 0.7], [3, 1, 2], [-1, 2, 3]], triangles[0]))
    R = (points[:, None, None] - triangles[None])[::-1, ::-1]
    normals, areas = normals[::-1], areas[::-1]
    edges = (np.roll(triangles, 1, -2) - np.roll(triangles, 2, -2))[::-1]
    for array in (R, normals, areas, edges):
        array.flags.writeable = False
    return R, normals, areas, edges


@pytest.mark.parametrize("supplied", [False, True])
@pytest.mark.parametrize("planar", [False, True])
def test_cached_triangle_formulas(geometry, supplied, planar):
    R, normals, areas, edges = geometry
    cache = TriangleIntegrals(
        R,
        normals,
        edges=edges if supplied else None,
        ta=areas if supplied else None,
        planar=planar,
    )
    references = {
        "gamma0": integrals.gamma0(R, edges=edges if supplied else None),
        "omega": integrals.omega(R),
        "x_distance": integrals.x_distance(R, normals),
        "normalized_x_distance": integrals.x_distance(R, normals, areas),
        "d_distance": integrals.d_distance(R, normals),
        "c_coeffs": integrals.c_coeffs(R, areas),
        "triangle_potential_uniform": integrals.triangle_potential_uniform(
            R, normals, planar, edges=edges if supplied else None
        ),
    }
    for name, expected in references.items():
        result = getattr(cache, name)
        np.testing.assert_allclose(result, expected, rtol=3e-12, atol=1e-12)
        assert getattr(cache, name) is result
    np.testing.assert_allclose(cache.ta, areas, rtol=1e-14)
    if supplied:
        assert cache.edges is edges
        assert cache.ta is areas
    with pytest.raises(FrozenInstanceError):
        cache.R = R.copy()


@pytest.mark.parametrize("reg,symmetrize", [(1e-13, True), (1e-9, False)])
def test_cached_edge_regularization(geometry, reg, symmetrize):
    R, normals, _, edges = geometry
    cache = TriangleIntegrals(R, normals, reg=reg, symmetrize=symmetrize, edges=edges)
    np.testing.assert_allclose(
        cache.gamma0,
        integrals.gamma0(R, reg=reg, symmetrize=symmetrize, edges=edges),
        rtol=1e-13,
        atol=1e-13,
    )


def test_distances_are_shared_and_planar_evaluation_is_lazy(geometry, monkeypatch):
    R, normals, _, edges = geometry
    calls = []
    original_norm = integrals_cached.norm

    def count_norm(value):
        calls.append(value)
        return original_norm(value)

    monkeypatch.setattr(integrals_cached, "norm", count_norm)
    cache = TriangleIntegrals(R, normals, edges=edges)
    _ = cache.gamma0
    _ = cache.omega
    _ = cache.triangle_potential_uniform
    assert sum(value is R for value in calls) == 1
    planar = TriangleIntegrals(R, normals, edges=edges, planar=True)
    _ = planar.triangle_potential_uniform
    assert "omega" not in planar.__dict__
    assert "determinant" not in planar.__dict__


def test_cached_dipole_and_approximate_potentials(geometry):
    R, normals, areas, edges = geometry
    # Keep approximate/dipole comparisons away from source vertices, where the
    # approximation has singular samples rather than an exact boundary formula.
    R = R[-3:]
    cache = TriangleIntegrals(R, normals, edges=edges, ta=areas)
    np.testing.assert_allclose(
        cache.triangle_potential_approx,
        integrals.triangle_potential_approx(R.mean(axis=-2), areas),
        rtol=1e-13,
    )
    np.testing.assert_allclose(
        cache.potential_dipoles,
        integrals.potential_dipoles(R, normals, areas),
        rtol=1e-13,
    )
    np.testing.assert_allclose(
        cache.triangle_potential_dipole_linear,
        integrals.triangle_potential_dipole_linear(R, normals, areas),
        rtol=1e-12,
        atol=1e-13,
    )


def test_cached_single_triangle(geometry):
    R, normals, areas, edges = geometry
    cache = TriangleIntegrals(R[:, 0], normals[0], edges=edges[0], ta=areas[0])
    expected = integrals.triangle_potential_uniform(R[:, :1], normals[:1])[:, 0]
    np.testing.assert_allclose(
        cache.triangle_potential_uniform, expected, rtol=1e-12, atol=1e-13
    )
    np.testing.assert_allclose(
        cache.potential_dipoles,
        integrals.potential_dipoles(R[:, :1], normals[:1], areas[:1])[:, 0],
        rtol=1e-13,
    )


def test_cached_vertex_dipoles():
    rng = np.random.default_rng(12)
    R = rng.normal(size=(5, 7, 3))
    normals, areas = rng.normal(size=(7, 3)), rng.random(7)
    cache = TriangleIntegrals(R, normals, ta=areas)
    np.testing.assert_allclose(
        cache.potential_vertex_dipoles,
        integrals.potential_vertex_dipoles(R, normals, areas),
        rtol=1e-14,
    )
