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


def test_mesh_projection_linear_function():
    import trimesh

    mesh = trimesh.Trimesh(
        vertices=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        faces=[[0, 1, 2]],
        process=False,
    )
    projection = utils.MeshProjection(mesh, 3)
    constant = projection.hatfunc_innerproducts(lambda points: np.ones(len(points)))
    linear = projection.hatfunc_innerproducts(lambda points: points[:, 0])
    np.testing.assert_allclose(constant, np.full(3, mesh.area / 3))
    np.testing.assert_allclose(linear, mesh.area * np.array([1.0, 2.0, 1.0]) / 12)
