from bfieldtools import quadratures
from bfieldtools import utils

from .test_line_conductor import _fake_line_conductor


import pytest
import numpy as np


@pytest.mark.parametrize(
    "rule", quadratures.dunavant_rules + quadratures.gauss_legendre_rules
)
def test_rules(rule):
    assert np.allclose(rule.weights.sum(), 1)
    assert np.allclose(rule.points.sum(axis=1), 1)


@pytest.mark.parametrize(
    "method", [("dunavant", 0), ("dunavant", 1), ("dunavant", 2), ("dunavant", 3)]
)
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


def test_dunavant3():
    p = np.array(
        [
            [1 - 2 * 0.445948490915965, 0.445948490915965, 0.445948490915965],
            [0.445948490915965, 1 - 2 * 0.445948490915965, 0.445948490915965],
            [0.445948490915965, 0.445948490915965, 1 - 2 * 0.445948490915965],
            [1 - 2 * 0.091576213509771, 0.091576213509771, 0.091576213509771],
            [0.091576213509771, 1 - 2 * 0.091576213509771, 0.091576213509771],
            [0.091576213509771, 0.091576213509771, 1 - 2 * 0.091576213509771],
        ]
    )
    pq = quadratures.get_triangle_rule(("dunavant", 3)).points
    assert np.allclose(pq, p)


def test_line_quad_points():
    lp = _fake_line_conductor()

    w, qp = quadratures.get_line_quad_points(lp.vertices[lp.entities[0].points])
    w0, qp0 = quadratures.get_line_quad_points(
        lp.vertices[lp.entities[0].points], "midpoint"
    )
    w1, qp1 = quadratures.get_line_quad_points(
        lp.vertices[lp.entities[0].points], ("gauss_legendre", 0)
    )

    assert np.allclose(w0, w1)
    assert np.allclose(qp0, qp1)

    w, qp = quadratures.get_line_quad_points(
        lp.vertices[lp.entities[0].points], "trapezoid"
    )
    w, qp = quadratures.get_line_quad_points(
        lp.vertices[lp.entities[0].points], "simpson"
    )
    w, qp = quadratures.get_line_quad_points(
        lp.vertices[lp.entities[0].points], ("gauss_legendre", 2)
    )


def test_spherical_quad_points():
    qp = quadratures.get_spherical_quad_points()
    assert np.allclose(sum(qp.weights), 1)
    for coord in qp.points.T:
        assert np.allclose(coord.sum(), 0, atol=1e-7)
