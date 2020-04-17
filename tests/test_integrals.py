from bfieldtools import integrals
from bfieldtools.mesh_calculus import mass_matrix
from bfieldtools.utils import load_example_mesh

import pytest

import numpy as np
from numpy.testing import assert_allclose
from scipy.sparse import csc_matrix


def test_gamma0():
    """
    No specific test done for this
    """

    pass


def test_omega():
    """
    Tests solid angle (whether sphere sums up to 4*pi)

    """

    mesh = load_example_mesh("unit_sphere")

    points = np.array([[0, 0, 0], [0.1, 0, 0.2], [-0.1, 0.1, -0.2]])

    R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]

    solid_angle = integrals.omega(R)

    assert_allclose(np.sum(solid_angle, axis=1), 4 * np.pi)


def test_x_distance():

    import trimesh

    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0]])

    faces = np.array([[2, 1, 0], [0, 3, 2]])

    mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

    points = np.array([[0, 1, 0]])

    R = points[:, None, None, :] - mesh.vertices[faces][None, :, :, :]

    x = integrals.x_distance(R, mesh.face_normals, mesh.area_faces)

    assert_allclose(x, np.array([[[0.0, 1.0, 0.0], [2.0, -1.0, -1.0]]]))


def test_d_distance():

    mesh = load_example_mesh("unit_disc")

    for z_offset in [1, 0, -1]:
        points = mesh.vertices + np.array([[0, 0, z_offset]])

        R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]

        assert_allclose(integrals.d_distance(R, mesh.face_normals), z_offset)


def test_dipole_potential_approximation():
    """
    Test whether approximations and exact solution are close enough

    """

    mesh = load_example_mesh("unit_disc")

    points = np.array([[0.1, 0, 0.2], [-0.1, 0.1, -0.2]]) * 20

    R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]
    R_v = points[:, None, :] - mesh.vertices[None, :, :]

    vertex_areas = mass_matrix(mesh, lumped=True).diagonal()

    approx_pot_v = integrals.potential_vertex_dipoles(
        R_v, mesh.vertex_normals, vertex_areas
    )

    approx_pot_f = integrals.potential_dipoles(R, mesh.face_normals, mesh.area_faces)
    exact_pot_f = integrals.triangle_potential_dipole_linear(
        R, mesh.face_normals, mesh.area_faces
    )

    # Map faces -> vertices
    Nf = len(mesh.faces)
    Nv = len(mesh.vertices)
    M0 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 0])), (Nf, Nv))
    M1 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 1])), (Nf, Nv))
    M2 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 2])), (Nf, Nv))
    exact_pot_v = (
        exact_pot_f[:, :, 0] @ M0
        + exact_pot_f[:, :, 1] @ M1
        + exact_pot_f[:, :, 2] @ M2
    )
    approx_pot_fv = (
        approx_pot_f[:, :, 0] @ M0
        + approx_pot_f[:, :, 1] @ M1
        + approx_pot_f[:, :, 2] @ M2
    )

    assert_allclose(approx_pot_v, exact_pot_v, rtol=5e-2)
    assert_allclose(approx_pot_fv, exact_pot_v, rtol=1e-3)


def test_triangle_potential_approximation():
    """
    test whether 1/r triangle potential approximations give close enough
    to exact (with/without planar assumption)
    """

    mesh = load_example_mesh("unit_disc")

    points = np.array([[0, 0, 1], [-0.1, 0.1, -5]])

    R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]
    Rcenters = points[:, None, :] - mesh.triangles_center[None, :, :]

    exact = integrals.triangle_potential_uniform(R, mesh.face_normals)
    approx = integrals.triangle_potential_approx(Rcenters, mesh.area_faces)

    assert_allclose(approx, exact, rtol=5e-3)

    points = np.array([[0, 0, 0], [-0.1, 0.1, 0]])
    R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]

    exact = integrals.triangle_potential_uniform(R, mesh.face_normals)
    exact_planar = integrals.triangle_potential_uniform(
        R, mesh.face_normals, planar=True
    )

    assert_allclose(exact_planar, exact)
