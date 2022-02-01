# -*- coding: utf-8 -*-
"""
Demonstration of contour manipulation using a Fourier basis
"""

import sys

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

    x_new = np.linspace(0, 2 * np.pi, 100)
    polys_new = []
    for polyline in polys:
        W = fit_fourier_basis(polyline, fourier_basis_funcs, freqs, 1e-4)

        polyline_new = (W.T @ fourier_basis_funcs(x_new)).T
        polys_new.append(polyline_new)
        plt.plot(polyline_new[:, 0], polyline_new[:, 1])

    for p in polys:
        plt.plot(p[:, 0], p[:, 1], "--*")

    plt.axis("equal")


def half_sigmoid(x, offset=0, scale=3):
    x = (x - 1) / (1 - offset)
    val = 1 / (1 + np.exp(-x * scale))
    val0 = 1 / (1 + np.exp(scale / (1 - offset)))
    val1 = 1 / 2
    return (val - val0) / (val1 - val0) * 0.5


def interp3(a, b, c, t, alpha=1, beta=1, tuning=0.75):
    tt1 = half_sigmoid(1 - t, np.sqrt(1 - alpha ** 2) * tuning)
    tt2 = half_sigmoid(t, np.sqrt(1 - beta ** 2) * tuning)

    f1 = a * tt1 + b * (1 - tt1)
    f2 = b * (1 - tt2) + c * tt2

    return (1 - t) * f1 + t * f2


def test_contour_merge(tuning=0.75):
    mesh, scalars = setup_contour_input()
    N = 10
    polys, vals = contour.scalar_contour(
        mesh, scalars, N_contours=N, return_values=True
    )
    fourier_basis_funcs, freqs = get_fourier_basis(4)

    x_new = np.linspace(0, 2 * np.pi, 500)
    Ws = []
    for polyline in polys:
        # polyline = polys[-1]
        W = fit_fourier_basis(polyline, fourier_basis_funcs, freqs, 1e-4)
        Ws.append(W)

    Ws_prev = Ws[:1] + Ws[:-1]
    Ws_next = Ws[1:] + Ws[-1:]
    for Wcurr, Wprev, Wnext in zip(Ws, Ws_prev, Ws_next):
        t_all = x_new / (2 * np.pi)
        polyline_new = np.zeros((len(x_new), 3))
        for ii, (f, t) in enumerate(zip(fourier_basis_funcs(x_new).T, t_all)):
            alpha = np.sum(Wcurr * Wprev) / np.sqrt(
                (np.sum(Wcurr ** 2) * np.sum(Wprev ** 2))
            )
            beta = np.sum(Wcurr * Wnext) / np.sqrt(
                (np.sum(Wcurr ** 2) * np.sum(Wnext ** 2))
            )
            W = interp3(Wprev, Wcurr, Wnext, t, alpha, beta, tuning=tuning)
            polyline_new[ii] = W.T @ f
        plt.plot(polyline_new[:, 0], polyline_new[:, 1])
    plt.axis("equal")


def test_contour_merge_fft():
    """
    Contour merge with global smoothing using FFT
    """
    mesh, scalars = setup_contour_input()
    N = 10
    polys, vals = contour.scalar_contour(
        mesh, scalars, N_contours=N, return_values=True
    )
    fourier_basis_funcs, freqs = get_fourier_basis(4)
    x_new = np.linspace(0, 2 * np.pi, 100)
    polys_new = []
    for polyline in polys:
        W = fit_fourier_basis(polyline, fourier_basis_funcs, freqs, 1e-4)
        polyline_new = (W.T @ fourier_basis_funcs(x_new)).T
        polys_new.append(polyline_new)

    polys_arr = np.concatenate(polys_new)
    polys_arr2 = np.concatenate(polys_new[::-1])
    polys_arr_concat = np.concatenate((polys_arr, polys_arr2))

    # Filter concatenated coordinate functions using FFT
    polys_fft = np.fft.fftshift(np.fft.fft(polys_arr_concat, axis=0), axes=0)
    N = polys_fft.shape[0]
    from scipy.signal.windows import hann

    # Tune this to regulate low-pass filtering
    # Maybe a more easily-tunable window would be better
    # or any other filtering method
    w = hann(N // 10)
    n_pad = (N - len(w)) // 2
    w = np.pad(w, (n_pad, n_pad))

    polys_merged = np.fft.ifft(np.fft.ifftshift(w[:, None] * polys_fft, axes=0), axis=0)

    plt.plot(polys_merged[: N // 2, 0], polys_merged[: N // 2, 1])
    plt.plot(polys_merged[N // 2 :, 0], polys_merged[N // 2 :, 1], "--")


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


def zero_phase_func(polyline):
    r = np.linalg.norm(polyline, axis=-1)
    return np.argmin(r)


def contour_basis_matrix(polyline, basis_funcs, zero_phase_index=None):
    edges_lens = np.linalg.norm(np.roll(polyline, -1, 0) - polyline, axis=-1)
    x = np.cumsum(edges_lens) / np.sum(edges_lens)
    x = np.roll(x, 1, 0) * 2 * np.pi
    x[0] = 0
    if zero_phase_index is not None:
        x -= x[zero_phase_index(polyline)]

    return basis_funcs(x)


def fit_fourier_basis(polyline, basis_funcs, freqs, _lambda=0.01):
    M = contour_basis_matrix(polyline, basis_funcs, zero_phase_func)
    DD = np.diag(freqs ** 4)
    MM = M @ M.T
    m0 = np.linalg.eigvalsh(MM)[-1]

    scaling = _lambda * m0 * 4
    W = np.linalg.solve(MM + scaling * DD, M @ polyline)

    return W


#%% Plot examples


def plot_stream_function():
    from bfieldtools.viz import plot_data_on_vertices

    mesh, scalars = setup_contour_input()
    plot_data_on_vertices(mesh, scalars, ncolors=32)


def plot_sigmoids():
    offsets = np.linspace(-1, 0.95, 10)
    t = np.linspace(0, 2, 100)
    for offset in offsets:
        plt.plot(t, half_sigmoid(t, offset))
    plt.xlabel("$t$")
    plt.ylabel("$\sigma(t)$")
    plt.plot([1, 1], [0, 1], "k--")
    plt.legend([f"$t_0=${t0:.2f}" for t0 in offsets])


def plot_interpolation_functions():
    t = np.linspace(0, 1, 100)
    styles = ["-.", "--", "-"]
    for alpha, style in zip([0.3, 0.75, 0.9], styles):
        for t_start in [-1, 0, 1]:
            plt.plot(t + t_start, interp3(1, 0, 0, t, alpha, alpha, 1), style)
            plt.plot(t + t_start, interp3(0, 1, 0, t, alpha, alpha, 1), style)
            plt.plot(t + t_start, interp3(0, 0, 1, t, alpha, alpha, 1), style)
            plt.gca().set_prop_cycle(None)
    plt.xlabel("$t$")
    plt.plot([1, 1], [0, 1], "k--")
    plt.plot([0, 0], [0, 1], "k--")
    plt.legend(("$h_{i-1}$", "$h_{i}$", "$h_{i+1}$"))


if __name__ == "__main__":
    test_contour_basis()
    # plt.figure()
    # test_contour_merge()
    plt.figure()
    test_contour_merge_fft()
    # plt.figure()
    # plot_sigmoids()
    # plt.figure()
    # plot_interpolation_functions()
    # plot_stream_function()
