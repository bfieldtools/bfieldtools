from bfieldtools import utils

import numpy as np

import pytest

from .test_line_conductor import _fake_line_conductor


def test_mesh_utils():

    mesh1 = utils.load_example_mesh("unit_disc")
    mesh2 = utils.load_example_mesh("10x10_plane")

    mesh2.vertices += np.array([[0, 0, 20]])

    mesh3 = utils.combine_meshes((mesh1, mesh2))

    assert len(mesh3.vertices) == len(mesh1.vertices) + len(mesh2.vertices)


@pytest.mark.xfail
def test_load_mesh_fail():
    mesh = utils.load_example_mesh("this mesh does not exist")


def test_quad_points():
    mesh = utils.load_example_mesh("unit_disc")

    w, qp = utils.get_quad_points(mesh.vertices, mesh.faces)
    w, qp = utils.get_quad_points(mesh.vertices, mesh.faces, "dunavant_02")
    w, qp = utils.get_quad_points(mesh.vertices, mesh.faces, method="lether", index=3)


@pytest.mark.xfail
def test_quad_points_fail1():
    mesh = utils.load_example_mesh("unit_disc")
    w, qp = utils.get_quad_points(mesh.vertices, mesh.faces, method="lether")


@pytest.mark.xfail
def test_quad_points_fail2():
    mesh = utils.load_example_mesh("unit_disc")
    w, qp = utils.get_quad_points(mesh.vertices, mesh.faces, method="dunavant")


def test_line_quad_points():
    lp = _fake_line_conductor()

    w, qp = utils.get_line_quad_points(lp.vertices[lp.entities[0].points])
    w, qp = utils.get_line_quad_points(lp.vertices[lp.entities[0].points], "midpoint")
    w, qp = utils.get_line_quad_points(lp.vertices[lp.entities[0].points], "fejer_1", 1)


@pytest.mark.xfail
def test_line_quad_points_fail1():
    lp = _fake_line_conductor()

    w, qp = utils.get_line_quad_points(lp.vertices[lp.entities[0].points], "fejer_1")


@pytest.mark.xfail
def test_line_quad_points_fail2():
    lp = _fake_line_conductor()

    w, qp = utils.get_line_quad_points(
        lp.vertices[lp.entities[0].points], method="dunavant"
    )


def test_cylinder_points():

    cp = utils.cylinder_points()
    cp = utils.cylinder_points(alpha=90)
    cp = utils.cylinder_points(orientation=np.array([0, 1, 0]))
