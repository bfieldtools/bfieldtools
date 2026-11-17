"""Independent kernel checks and physical convergence of constant-face BEM."""

import numpy as np
import pytest
import trimesh
from numpy.testing import assert_allclose
from scipy.integrate import quad_vec

from bfieldtools.bem import (
    MagneticBoundary,
    adjoint_double_layer,
    single_layer_gradient,
    single_layer_potential,
)


@pytest.fixture
def triangle():
    return trimesh.Trimesh(
        vertices=[[0, 0, 0], [1.1, 0.1, 0], [-0.1, 0.9, 0]],
        faces=[[0, 1, 2]],
        process=False,
    )


def test_integrals_against_quadrature(triangle):
    """Independent integration in barycentric coordinates, both sides of a face."""
    points = np.array([[0.2, 0.3, 0.4], [1.8, -0.2, -0.6]])
    v0, v1, v2 = triangle.triangles[0]
    e1, e2 = v1 - v0, v2 - v0
    jacobian = np.linalg.norm(np.cross(e1, e2))

    def integrand(s, t):
        r = points - (v0 + s * e1 + t * e2)
        distance = np.linalg.norm(r, axis=1)
        return (
            jacobian
            / (4 * np.pi)
            * np.column_stack([1 / distance, -r / distance[:, None] ** 3])
        )

    reference, _ = quad_vec(
        lambda t: quad_vec(lambda s: integrand(s, t), 0, 1 - t)[0], 0, 1
    )
    assert_allclose(
        single_layer_potential(triangle, points)[:, 0], reference[:, 0], rtol=1e-10
    )
    assert_allclose(
        single_layer_gradient(triangle, points)[:, :, 0], reference[:, 1:], rtol=1e-10
    )


def test_gradient_and_scaling(triangle):
    points = np.array([[0.2, 0.3, 0.4], [1.8, -0.2, -0.6]])
    gradient = single_layer_gradient(triangle, points)
    transform = trimesh.transformations.rotation_matrix(0.7, [1, 2, 3])
    transform[:3, 3] = [0.2, -0.4, 1.3]
    rotated = triangle.copy()
    rotated.apply_transform(transform)
    transformed_points = points @ transform[:3, :3].T + transform[:3, 3]
    assert_allclose(
        single_layer_potential(rotated, transformed_points),
        single_layer_potential(triangle, points),
        atol=1e-12,
    )
    assert_allclose(
        single_layer_gradient(rotated, transformed_points),
        np.einsum("ab,pbf->paf", transform[:3, :3], gradient),
        atol=1e-12,
    )
    flipped = triangle.copy()
    flipped.invert()
    assert_allclose(single_layer_gradient(flipped, points), gradient, atol=1e-12)
    for a in range(3):
        step = np.eye(3)[a] * 1e-5
        difference = (
            single_layer_potential(triangle, points + step)
            - single_layer_potential(triangle, points - step)
        ) / 2e-5
        assert_allclose(difference, gradient[:, a], atol=1e-10)
    for scale in [1e-6, 1e6]:
        scaled = triangle.copy()
        scaled.apply_scale(scale)
        assert_allclose(
            single_layer_potential(scaled, points * scale),
            single_layer_potential(triangle, points) * scale,
            rtol=1e-10,
        )
        assert_allclose(
            single_layer_gradient(scaled, points * scale), gradient, atol=1e-12
        )


def test_jump_and_surface_rejection(triangle):
    center = triangle.triangles_center
    normal = triangle.face_normals[0]
    above = single_layer_gradient(triangle, center + 1e-7 * normal)[0, :, 0]
    below = single_layer_gradient(triangle, center - 1e-7 * normal)[0, :, 0]
    assert_allclose(below - above, normal, atol=1e-6)
    assert np.isfinite(single_layer_potential(triangle, center)).all()
    for point in [
        center,
        triangle.vertices[:1],
        triangle.vertices[:2].mean(axis=0)[None],
    ]:
        with pytest.raises(ValueError, match="off the source surface"):
            single_layer_gradient(triangle, point)


def test_chunking_and_empty_points():
    mesh = trimesh.creation.icosphere(subdivisions=1)
    points = np.array([[0.1, 0.2, 0.3], [2, 1, -0.2], [0, 0, 0]])
    for operator, shape in [
        (single_layer_potential, (0, len(mesh.faces))),
        (single_layer_gradient, (0, 3, len(mesh.faces))),
    ]:
        assert_allclose(
            operator(mesh, points, 1), operator(mesh, points, 5), atol=1e-14
        )
        assert operator(mesh, np.empty((0, 3))).shape == shape
    boundary = MagneticBoundary(mesh)
    assert_allclose(
        adjoint_double_layer(boundary, 1), adjoint_double_layer(boundary, 5), atol=1e-14
    )


def test_boundary_shell_and_material_mapping():
    outer = trimesh.creation.icosphere(subdivisions=1)
    inner = trimesh.creation.icosphere(subdivisions=1, radius=0.6)
    inner.invert()
    shell = trimesh.util.concatenate([outer, inner])
    other = trimesh.creation.icosphere(subdivisions=0, radius=0.2)
    other.apply_translation([3, 0, 0])
    boundary = MagneticBoundary([shell, other], mu_r=[10, 20])
    assert_allclose(boundary.mesh.triangles[: len(shell.faces)], shell.triangles)
    assert_allclose(boundary.mesh.face_normals[: len(shell.faces)], shell.face_normals)
    assert_allclose(boundary.face_mu_r[: len(shell.faces)], 10)
    assert_allclose(boundary.beta[-len(other.faces) :], 19 / 21)
    assert len(np.unique(boundary.component_ids)) == 3
    assert_allclose(boundary.body_ids[-len(other.faces) :], 1)
    # A separate solid in the vacuum cavity is also supported.
    island = trimesh.creation.icosphere(subdivisions=0, radius=0.2)
    MagneticBoundary([shell, island], mu_r=[10, 20])
    inner.invert()
    with pytest.raises(ValueError, match="Normals"):
        MagneticBoundary(trimesh.util.concatenate([outer, inner]))


@pytest.mark.parametrize("mu_r", [0, -1, np.inf, np.nan, [2, 3]])
def test_invalid_material(mu_r):
    with pytest.raises(ValueError, match="mu_r"):
        MagneticBoundary(trimesh.creation.box(), mu_r)


def test_invalid_geometry():
    box = trimesh.creation.box()
    with pytest.raises(ValueError, match="watertight"):
        MagneticBoundary(trimesh.Trimesh(box.vertices, box.faces[:-1], process=False))
    inward = box.copy()
    inward.invert()
    with pytest.raises(ValueError, match="Normals"):
        MagneticBoundary(inward)
    for displacement in [[0.5, 0.1, 0.2], [1, 0, 0], [0, 0, 0]]:
        other = box.copy()
        other.apply_translation(displacement)
        with pytest.raises(ValueError, match="intersect or touch"):
            MagneticBoundary([box, other])
    with pytest.raises(ValueError, match="Normals"):
        MagneticBoundary([box, trimesh.creation.box(extents=[0.5] * 3)])
    other = box.copy()
    other.apply_translation([3, 0, 0])
    with pytest.raises(ValueError, match="one connected metal body"):
        MagneticBoundary(trimesh.util.concatenate([box, other]))


def test_sphere_convergence():
    errors = []
    for refinement in [1, 2, 3]:
        mesh = trimesh.creation.icosphere(subdivisions=refinement)
        boundary = MagneticBoundary(mesh, mu_r=10)
        k = adjoint_double_layer(boundary)
        assert_allclose(np.diag(k), 0, atol=0)
        g = mesh.face_normals[:, 2]
        q = np.linalg.solve(
            np.eye(len(g)) + 2 * boundary.beta[:, None] * k, 2 * boundary.beta * g
        )
        h = np.array([0, 0, 1.0]) - single_layer_gradient(mesh, [[0, 0, 0]]) @ q
        errors.append(abs(h[0, 2] - 3 / 12))
        assert_allclose(np.sum(mesh.area_faces * q), 0, atol=1e-12)
        # No contrast must give zero response, using the same geometry operator.
        beta = MagneticBoundary(mesh, mu_r=1).beta
        q0 = np.linalg.solve(np.eye(len(g)) + 2 * beta[:, None] * k, 2 * beta * g)
        assert_allclose(q0, 0, atol=0)
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < 0.04  # Smooth-sphere error includes faceting and collocation.


def test_close_shell_blocks():
    """Analytic cross-surface kernels retain the jump across a very thin gap."""
    outer = trimesh.creation.icosphere(subdivisions=1)
    inner = outer.copy()
    inner.apply_scale(1 - 1e-5)
    inner.invert()
    boundary = MagneticBoundary(trimesh.util.concatenate([outer, inner]))
    k = adjoint_double_layer(boundary)
    n = len(outer.faces)
    assert np.isfinite(k).all()
    # Corresponding inner/outer triangles become parallel coincident panels.
    assert_allclose(np.diag(k[:n, n:]), -0.5, atol=1e-4)
    assert_allclose(np.diag(k[n:, :n]), -0.5, atol=1e-4)


def test_shell_cavity_convergence():
    """Two-surface solve converges to the analytic spherical-shell cavity field."""
    mu_r, inner_radius = 10.0, 0.6
    exact = (
        9
        * mu_r
        / ((2 * mu_r + 1) * (mu_r + 2) - 2 * inner_radius**3 * (mu_r - 1) ** 2)
    )
    errors = []
    for refinement in [1, 2]:
        outer = trimesh.creation.icosphere(subdivisions=refinement)
        inner = outer.copy()
        inner.apply_scale(inner_radius)
        inner.invert()
        boundary = MagneticBoundary(trimesh.util.concatenate([outer, inner]), mu_r)
        k = adjoint_double_layer(boundary)
        g = boundary.mesh.face_normals[:, 2]
        q = np.linalg.solve(
            np.eye(len(g)) + 2 * boundary.beta[:, None] * k, 2 * boundary.beta * g
        )
        h = (
            np.array([0.0, 0, 1])
            - single_layer_gradient(boundary.mesh, [[0, 0, 0]]) @ q
        )
        errors.append(abs(h[0, 2] - exact))
    assert errors[1] < errors[0]
    assert errors[1] < 0.12
