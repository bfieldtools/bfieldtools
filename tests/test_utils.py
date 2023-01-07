from bfieldtools import utils

import numpy as np

import pytest

def test_mesh_utils():

    mesh1 = utils.load_example_mesh("unit_disc")
    mesh2 = utils.load_example_mesh("10x10_plane")

    mesh2.vertices += np.array([[0, 0, 20]])

    mesh3 = utils.combine_meshes((mesh1, mesh2))

    assert len(mesh3.vertices) == len(mesh1.vertices) + len(mesh2.vertices)


@pytest.mark.xfail
def test_load_mesh_fail():
    mesh = utils.load_example_mesh("this mesh does not exist")


def test_cylinder_points():

    cp = utils.cylinder_points()
    cp = utils.cylinder_points(alpha=90)
    cp = utils.cylinder_points(orientation=np.array([0, 1, 0]))
