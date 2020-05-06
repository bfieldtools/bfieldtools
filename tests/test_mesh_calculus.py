from bfieldtools import mesh_calculus
from bfieldtools.utils import load_example_mesh

import pytest

import numpy as np
from numpy.testing import assert_allclose


def test_laplacian():

    mesh = load_example_mesh("unit_disc")

    L = mesh_calculus.laplacian_matrix(mesh)


def test_mass():

    mesh = load_example_mesh("unit_disc")

    M = mesh_calculus.mass_matrix(mesh, lumped=True)
    Ml = mesh_calculus.mass_matrix(mesh, lumped=False)


def test_others():

    mesh = load_example_mesh("unit_disc")

    vals = np.ones((len(mesh.vertices),)) - np.linalg.norm(mesh.vertices, axis=1)

    G = mesh_calculus.gradient(vals, mesh, rotated=False)
    Gr = mesh_calculus.gradient(vals, mesh, rotated=True)

    D = mesh_calculus.divergence(Gr.T, mesh)
    C = mesh_calculus.curl(Gr.T, mesh)
