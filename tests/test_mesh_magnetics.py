from bfieldtools import mesh_magnetics
from bfieldtools.utils import load_example_mesh

import pytest

import numpy as np
from numpy.testing import assert_allclose


def test_magnetic_field_coupling():

    mesh = load_example_mesh("unit_disc")

    points = mesh.vertices + np.array([0, 0, 20])

    B_1_chunk = mesh_magnetics.magnetic_field_coupling(mesh, points, Nchunks=1)
    B_2_chunk = mesh_magnetics.magnetic_field_coupling(mesh, points, Nchunks=2)

    assert_allclose(B_1_chunk, B_2_chunk, atol=1e-16)

    B_q_1 = mesh_magnetics.magnetic_field_coupling(mesh, points, quad_degree=1)
    B_q_2 = mesh_magnetics.magnetic_field_coupling(mesh, points, quad_degree=2)
    B_q_3 = mesh_magnetics.magnetic_field_coupling(mesh, points, quad_degree=3)

    assert_allclose(B_q_1, B_q_2, atol=1e-3 * np.mean(np.abs(B_q_2)))
    assert_allclose(B_q_2, B_q_3, atol=1e-3 * np.mean(np.abs(B_q_2)))

    B_a = mesh_magnetics.magnetic_field_coupling(mesh, points, analytic=True)

    assert_allclose(B_a, B_q_1, atol=1e-3 * np.mean(np.abs(B_a)))
    assert_allclose(B_a, B_q_2, atol=1e-3 * np.mean(np.abs(B_a)))
    assert_allclose(B_a, B_q_3, atol=1e-3 * np.mean(np.abs(B_a)))


def test_magnetic_scalar_potential_coupling():

    mesh = load_example_mesh("unit_disc")

    points = mesh.vertices + np.array([0, 0, 20])

    U_1_chunk = mesh_magnetics.scalar_potential_coupling(mesh, points, Nchunks=1)
    U_2_chunk = mesh_magnetics.scalar_potential_coupling(mesh, points, Nchunks=2)

    assert_allclose(U_1_chunk, U_2_chunk, atol=1e-16)

    U_exact = mesh_magnetics.scalar_potential_coupling(mesh, points, approx_far=False)
    U_approx_far = mesh_magnetics.scalar_potential_coupling(
        mesh, points, approx_far=True, margin=0.5
    )

    assert_allclose(U_exact, U_approx_far, atol=2e-3 * np.mean(np.abs(U_exact)))


def test_magnetic_vector_potential_coupling():

    mesh = load_example_mesh("unit_disc")

    points = mesh.vertices + np.array([0, 0, 20])

    A_1_chunk = mesh_magnetics.vector_potential_coupling(mesh, points, Nchunks=1)
    A_2_chunk = mesh_magnetics.vector_potential_coupling(mesh, points, Nchunks=2)

    assert_allclose(A_1_chunk, A_2_chunk, atol=1e-16)

    A_exact = mesh_magnetics.vector_potential_coupling(mesh, points, approx_far=False)
    A_approx_far = mesh_magnetics.vector_potential_coupling(
        mesh, points, approx_far=True, margin=0.5
    )
    A_approx_far_clusters = mesh_magnetics.vector_potential_coupling(
        mesh, points, approx_far=True, margin=0.5, chunk_clusters=True
    )

    assert_allclose(A_exact, A_approx_far, atol=2e-3 * np.mean(np.abs(A_exact)))
    assert_allclose(
        A_exact, A_approx_far_clusters, atol=2e-3 * np.mean(np.abs(A_exact))
    )
