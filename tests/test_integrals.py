# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 11:32:38 2020

@author: Rasmus Zetter
"""


from bfieldtools.integrals import (
    omega,
    d_distance,
    potential_vertex_dipoles,
    triangle_potential_dipole_linear,
    potential_dipoles,
)
from bfieldtools.mesh_calculus import mass_matrix
from bfieldtools.utils import load_example_mesh


import pytest

import numpy as np
from numpy.testing import assert_allclose
from scipy.sparse import csc_matrix


def test_solid_angle():
    mesh = load_example_mesh("unit_sphere")

    points = np.array([[0, 0, 0], [0.1, 0, 0.2], [-0.1, 0.1, -0.2]])

    R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]

    solid_angle = omega(R)

    assert_allclose(np.sum(solid_angle, axis=1), 4 * np.pi)


def test_d_distance():

    mesh = load_example_mesh("unit_disc")

    for z_offset in [1, 0, -1]:
        points = mesh.vertices + np.array([[0, 0, z_offset]])

        R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]

        assert_allclose(d_distance(R, mesh.face_normals), z_offset)


def test_dipole_potential_approximation():

    mesh = load_example_mesh("unit_disc")

    points = np.array([[0.1, 0, 0.2], [-0.1, 0.1, -0.2]]) * 20

    R = points[:, None, None, :] - mesh.vertices[mesh.faces][None, :, :, :]
    R_v = points[:, None, :] - mesh.vertices[None, :, :]

    vertex_areas = mass_matrix(mesh, lumped=True).diagonal()

    approx_pot_v = potential_vertex_dipoles(R_v, mesh.vertex_normals, vertex_areas)

    approx_pot_f = potential_dipoles(R, mesh.face_normals, mesh.area_faces)
    exact_pot_f = triangle_potential_dipole_linear(
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
