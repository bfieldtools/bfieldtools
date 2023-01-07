"""

Quadrature rules for numerical integration

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
)


def get_rule(method):
    name, degree = method
    rules = {"dunavant": dunavant_rules}
    return rules[name][degree]


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

    rule = get_rule(method)
    qp = np.einsum("ij,kjl->kil", rule.points, verts[tris])

    if return_ref_coords:
        return rule.weights, qp, rule.points

    return rule.weights, qp
