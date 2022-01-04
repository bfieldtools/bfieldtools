# -*- coding: utf-8 -*-
"""
Created on Sun May  9 16:20:16 2021

@author: antti
"""

import sys

sys.path.insert(0, "C:/Users/antti/bfieldtools")
from bfieldtools import contour
import numpy as np
import matplotlib.pyplot as plt


def setup_contour_input():
    """ Load example mesh and create scalars data
    """
    from bfieldtools.utils import load_example_mesh

    mesh = load_example_mesh("unit_disc")

    r = np.linalg.norm(mesh.vertices, axis=1)
    scalars = (1 - r) ** 2
    scalars *= mesh.vertices[:, 0]

    return mesh, scalars


def test_contour_basis():
    mesh, scalars = setup_contour_input()
    N = 10
    polys, vals = contour.scalar_contour(
        mesh, scalars, N_contours=N, return_values=True
    )
    fourier_basis_funcs, freqs = get_fourier_basis(4)

    for polyline in polys:
        # polyline = polys[-1]
        W = fit_fourier_basis(polyline, fourier_basis_funcs, freqs, 1e-4)

        x_new = np.linspace(0, 2 * np.pi, 100)
        polyline_new = (W.T @ fourier_basis_funcs(x_new)).T
        plt.plot(polyline_new[:, 0], polyline_new[:, 1])

    for p in polys:
        plt.plot(p[:, 0], p[:, 1], "--*")

    plt.axis("equal")


def get_fourier_basis(Nfreqs):
    freqs = np.zeros(2 * Nfreqs + 1)
    for k in range(1, Nfreqs + 1):
        freqs[2 * k - 1] = k
        freqs[2 * k] = k

    def basis_funcs(x):
        M = np.zeros((2 * Nfreqs + 1, len(x)))
        M[0] = 1
        for k in range(1, Nfreqs + 1):
            M[2 * k - 1] = np.cos(k * x)
            M[2 * k] = np.sin(k * x)
        return M

    return basis_funcs, freqs


def contour_basis_matrix(polyline, basis_funcs):
    edges_lens = np.linalg.norm(np.roll(polyline, -1, 0) - polyline, axis=-1)
    x = np.cumsum(edges_lens) / np.sum(edges_lens)
    x = np.roll(x, 1, 0) * 2 * np.pi
    x[0] = 0

    return basis_funcs(x)


def fit_fourier_basis(polyline, basis_funcs, freqs, _lambda=0.01):
    M = contour_basis_matrix(polyline, basis_funcs)
    DD = np.diag(freqs ** 4)
    MM = M @ M.T
    m0 = np.linalg.eigvalsh(MM)[-1]

    scaling = _lambda * m0 * 4
    W = np.linalg.solve(MM + scaling * DD, M @ polyline)

    return W


if __name__ == "__main__":
    test_contour_basis()
