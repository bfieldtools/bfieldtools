"""Geometry and optional-backend integration checks for volume meshing."""

import sys

import numpy as np
import pytest
import trimesh

from bfieldtools.utils import load_example_mesh
from bfieldtools.volume import TetraMesh, tetrahedralize


def test_tetrahedron_geometry():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1.0]])
    cells = np.array([[1, 0, 2, 3]])  # Deliberately inverted.
    mesh = TetraMesh(vertices, cells)
    np.testing.assert_allclose(mesh.volumes, [1 / 6])
    np.testing.assert_allclose(mesh.cell_centers, [[0.25, 0.25, 0.25]])
    assert cells[0, 0] == 1  # Caller data remains unchanged.
    triangles = mesh.vertices[mesh.boundary_faces]
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    assert np.all(np.einsum("ij,ij->i", normals, triangles.mean(axis=1) - 0.25) > 0)
    assert not mesh.vertices.flags.writeable
    grid = mesh.to_pyvista()
    assert grid.n_cells == 1
    np.testing.assert_allclose(grid.compute_cell_sizes()["Volume"], mesh.volumes)


def test_regular_tetrahedron_quality():
    mesh = TetraMesh([[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]], [[0, 1, 2, 3]])
    np.testing.assert_allclose(mesh.cell_quality, 1)


def test_sliver_rejected_before_integration():
    vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0.3, 0.3, 1e-5]]
    with pytest.raises(ValueError, match="tetrahedron 0.*quality.*min_quality"):
        TetraMesh(vertices, [[0, 1, 2, 3]])
    mesh = TetraMesh(vertices, [[0, 1, 2, 3]], min_quality=0)
    assert 0 < mesh.cell_quality[0] < 0.01
    with pytest.raises(ValueError, match="quality"):
        TetraMesh(vertices, [[0, 1, 2, 3]], min_quality=0.5)


@pytest.mark.parametrize("minimum", [-1, 1.1, np.nan, np.inf])
def test_invalid_minimum_quality(minimum):
    with pytest.raises(ValueError, match="min_quality"):
        TetraMesh(
            np.vstack((np.zeros(3), np.eye(3))), [[0, 1, 2, 3]], min_quality=minimum
        )


@pytest.mark.parametrize(
    "cells",
    [
        [[0, 1, 2, 2]],
        [[0, 1, 2, 4]],
        [[0.0, 1.0, 2.0, 3.0]],
        [[0, 1, 2, 3], [3, 2, 1, 0]],
    ],
)
def test_invalid_connectivity(cells):
    with pytest.raises(ValueError):
        TetraMesh(np.eye(4, 3), cells)


def test_degenerate_cell():
    with pytest.raises(ValueError, match="degenerate"):
        TetraMesh([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], [[0, 1, 2, 3]])


def test_shared_face():
    vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1]]
    mesh = TetraMesh(vertices, [[0, 1, 2, 3], [0, 1, 2, 4]])
    assert len(mesh.boundary_faces) == 6
    vertices[-1] = [0, 0, 2]  # Both cells now lie on the same side.
    with pytest.raises(ValueError, match="opposite"):
        TetraMesh(vertices, [[0, 1, 2, 3], [0, 1, 2, 4]])


def test_invalid_surface():
    sphere = trimesh.creation.icosphere(subdivisions=1)
    open_surface = trimesh.Trimesh(sphere.vertices, sphere.faces[:-1], process=False)
    with pytest.raises(ValueError, match="watertight"):
        tetrahedralize(open_surface)
    sphere.invert()
    with pytest.raises(ValueError, match="outward"):
        tetrahedralize(sphere)
    with pytest.raises(ValueError, match="single genus-zero"):
        tetrahedralize(
            trimesh.util.concatenate(
                [sphere, sphere.copy().apply_translation([4, 0, 0])]
            )
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_volume": 0},
        {"max_volume": np.nan},
        {"min_ratio": 1},
        {"min_dihedral": -1},
        {"min_dihedral": 60},
    ],
)
def test_invalid_options(kwargs):
    with pytest.raises(ValueError):
        tetrahedralize(trimesh.creation.icosphere(subdivisions=1), **kwargs)


def test_missing_backend(monkeypatch):
    monkeypatch.setitem(sys.modules, "tetgen", None)
    with pytest.raises(ImportError, match=r"bfieldtools\[meshing\]"):
        tetrahedralize(trimesh.creation.icosphere(subdivisions=1))


def test_sphere_tetgen():
    pytest.importorskip("tetgen")
    surface = load_example_mesh("unit_sphere")
    original_vertices = surface.vertices.copy()
    original_faces = surface.faces.copy()
    coarse = tetrahedralize(surface, max_volume=0.03)
    fine = tetrahedralize(surface, max_volume=0.005)
    for mesh in (coarse, fine):
        assert np.all(mesh.volumes > 0)
        np.testing.assert_allclose(mesh.volumes.sum(), surface.volume, rtol=1e-8)
        boundary = trimesh.Trimesh(mesh.vertices, mesh.boundary_faces, process=False)
        assert boundary.is_watertight and boundary.is_winding_consistent
        assert boundary.euler_number == 2
        assert len(boundary.split()) == 1
        np.testing.assert_allclose(boundary.volume, surface.volume, rtol=1e-8)
        np.testing.assert_allclose(boundary.area, surface.area, rtol=1e-8)
    assert len(fine.tetrahedra) > len(coarse.tetrahedra)
    assert fine.volumes.max() < coarse.volumes.max()
    np.testing.assert_array_equal(surface.vertices, original_vertices)
    np.testing.assert_array_equal(surface.faces, original_faces)
