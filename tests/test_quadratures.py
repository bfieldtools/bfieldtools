from bfieldtools import quadratures
from bfieldtools import utils

from .test_line_conductor import _fake_line_conductor


import pytest
import numpy as np


@pytest.mark.parametrize("rule", quadratures.dunavant_rules)
def test_rules(rule):
    assert np.allclose(rule.weights.sum(), 1)
    assert np.allclose(rule.points.sum(axis=1), 1)


@pytest.mark.parametrize("method", [("dunavant", 0), ("dunavant", 1), ("dunavant", 2)])
def test_test_integration(method):
    val = 5
    verts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    tris = np.array([[0, 1, 2]])

    def uniform(r):
        return val * np.ones(r.shape[:2])

    w, qp = quadratures.get_quad_points(verts, tris, method)
    assert np.allclose(uniform(qp) @ w, val)


def test_quad_points():
    mesh = utils.load_example_mesh("unit_disc")
    w, qp = quadratures.get_quad_points(mesh.vertices, mesh.faces, ("dunavant", 2))


def test_line_quad_points():
    lp = _fake_line_conductor()

    w, qp = quadratures.get_line_quad_points(lp.vertices[lp.entities[0].points])
    w, qp = quadratures.get_line_quad_points(lp.vertices[lp.entities[0].points], "midpoint")
    w, qp = quadratures.get_line_quad_points(lp.vertices[lp.entities[0].points], "trapezoid")
    w, qp = quadratures.get_line_quad_points(lp.vertices[lp.entities[0].points], "simpson")