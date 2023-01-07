"""

Quadrature rules for numerical integration

Refs:
D.A. Dunavant https://doi.org/10.1002/nme.1620210612

"""

import numpy as np
from dataclasses import dataclass


@dataclass
class QuadratureRule:
    points: np.ndarray
    weights: np.ndarray

    def __init__(self, points, weights):
        self.points = np.asarray(points)
        self.weights = np.asarray(weights)


dunavant_rules = (
    QuadratureRule([[1 / 3, 1 / 3, 1 / 3]], [1]),
    QuadratureRule(
        [[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]],
        [1 / 3] * 3,
    ),
    QuadratureRule(
        [
            [1 / 3, 1 / 3, 1 / 3],
            [3 / 5, 1 / 5, 1 / 5],
            [1 / 5, 3 / 5, 1 / 5],
            [1 / 5, 1 / 5, 3 / 5],
        ],
        [-0.5625] + [0.520833333333333] * 3,
    ),
    QuadratureRule(
        [
            [1 - 2 * 0.445948490915965, 0.445948490915965, 0.445948490915965],
            [0.445948490915965, 1 - 2 * 0.445948490915965, 0.445948490915965],
            [0.445948490915965, 0.445948490915965, 1 - 2 * 0.445948490915965],
            [1 - 2 * 0.091576213509771, 0.091576213509771, 0.091576213509771],
            [0.091576213509771, 1 - 2 * 0.091576213509771, 0.091576213509771],
            [0.091576213509771, 0.091576213509771, 1 - 2 * 0.091576213509771],
        ],
        [0.223381589678011] * 3 + [0.109951743655322] * 3,
    ),
)

line_rules = {
    "midpoint": QuadratureRule([[0.5, 0.5]], [1]),
    "trapezoid": QuadratureRule([[1, 0], [0, 1]], [0.5, 0.5]),
    "simpson": QuadratureRule([[1, 0], [0.5, 0.5], [0, 1]], [1 / 6, 4 / 6, 1 / 6]),
}


def get_triangle_rule(method):
    name, degree = method
    rules = {"dunavant": dunavant_rules}
    return rules[name][degree]


def get_line_rule(method):
    return line_rules[method]


def get_quad_points(verts, tris, method, return_ref_coords=False):
    """Get quad points and weights from quadrature rules implemented in
    quadpy  

    Parameters
    ----------
    verts: array-like [Nverts x 3]
    tris: array-like [Ntris x 3]
    method: tuple (method_name: str, degree: int)

    Returns
    -------
    w: array-like  (Nquad, )
        quadrature weights
    qp: array-like (Ntris, Nquad, xyz)
        quadrature points in each triangle

    """

    rule = get_triangle_rule(method)
    qp = np.einsum("ij,kjl->kil", rule.points, verts[tris])

    if return_ref_coords:
        return rule.weights, qp, rule.points

    return rule.weights, qp


def get_line_quad_points(line_vertices, method="midpoint", index=None):
    """Get quad points and weights from quadrature rules implemented in
    quadpy

    Parameters
    ----------
    line_vertices: array-like [Nverts x 3]
        Assumes vertices are connected according to index, last index connects to first
    method: str
        One of 'midpoint', 'trapezoid', or 'simpson'

    Returns
    -------
    w: array-like  (Nquad, )
        quadrature weights
    qp: array-like (Nverts, Nquad)
        quadrature points for each edge connecting the vertices

    """
    rule = get_line_rule(method)

    qp = np.zeros((len(line_vertices), len(rule.weights), 3))

    p0 = line_vertices
    p1 = np.roll(line_vertices, -1)

    B = np.array([p0, p1])
    qp = np.einsum("ij,jkl->kil", rule.points, B)

    return rule.weights, qp
