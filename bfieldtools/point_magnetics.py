# -*- coding: utf-8 -*-
"""
Fields for magnetic dipoles
"""

import numpy as np

from numba import jit


def cross(p, m):
    c = np.zeros_like(p)
    c[:, 0] = p[:, 1] * m[2] - p[:, 2] * m[1]
    c[:, 1] = p[:, 2] * m[0] - p[:, 0] * m[2]
    c[:, 2] = p[:, 0] * m[1] - p[:, 1] * m[0]

    return c


@jit
def magnetic_field_dipole_origin(points, moment):
    coeff = 1e-7  # mu0/(4*pi)
    m = moment
    p = points
    r = np.linalg.norm(points, axis=1)
    b = 3 * np.einsum("ij,ik,k->ij", p, p, m) / (r ** 2)[:, None] - m

    return coeff * b / (r ** 3)[:, None]


def vector_potential_dipole_origin(points, moment):
    coeff = 1e-7  # mu0/(4*pi)
    m = moment
    p = points
    r = np.linalg.norm(points, axis=1)
    return coeff * cross(p, m) / (r ** 3)[:, None]


@jit
def magnetic_field_dipoles(points, moments, points_dipoles):
    bb = np.zeros(moments.shape[:1] + points.shape, dtype=float)
    for ii, (p0, m) in enumerate(zip(points_dipoles, moments)):
        bb[ii] = magnetic_field_dipole_origin(points - p0, m)
    return bb


def vector_potential_dipoles(points, moments, points_dipoles):
    aa = np.zeros(moments.shape[:1] + points.shape, dtype=float)
    for ii, (p0, m) in enumerate(zip(points_dipoles, moments)):
        aa[ii] = vector_potential_dipole_origin(points - p0, m)
    return aa
