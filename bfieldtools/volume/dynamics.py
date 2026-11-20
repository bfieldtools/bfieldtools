"""External-loop coupling and modal evolution of conserved volume currents."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.constants import mu_0
from scipy.integrate import solve_ivp

from ..line_magnetics import vector_potential
from .basis import CellCurrentBasis
from .integrals import tetrahedron_quadrature, _positive_integer


def loop_coupling(
    basis: CellCurrentBasis,
    vertices: ArrayLike,
    *,
    quadrature_order: int = 4,
    cell_chunk_size: int = 32,
) -> NDArray[np.float64]:
    """Return m_a = integral j_a dot A_loop^(1) dV, in henries.

    vertices is a closed oriented polygon (first point repeated at the end),
    outside the conductor. Segments touching or crossing any tetrahedron are
    rejected. Exact straight-segment vector potentials are volume-integrated
    with positive tetrahedral quadrature. Refine the polygon separately from
    the volume quadrature when approximating a curved wire.
    """
    vertices = np.asarray(vertices)
    if np.iscomplexobj(vertices):
        raise ValueError("loop vertices must be real")
    vertices = np.asarray(vertices, dtype=float)
    if (
        vertices.ndim != 2
        or vertices.shape[1] != 3
        or len(vertices) < 4
        or not np.isfinite(vertices).all()
        or not np.array_equal(vertices[0], vertices[-1])
        or np.any(np.linalg.norm(np.diff(vertices, axis=0), axis=1) == 0)
    ):
        raise ValueError(
            "vertices must define a finite closed polygon with nonzero segments"
        )
    _positive_integer(cell_chunk_size, "cell_chunk_size")
    barycentric, weights = tetrahedron_quadrature(quadrature_order)
    mesh = basis.mesh
    cells = mesh.vertices[mesh.tetrahedra]
    transforms = np.linalg.inv(np.swapaxes(cells[:, 1:] - cells[:, :1], 1, 2))
    # Clip each segment to all four barycentric half-spaces. This also catches
    # wires passing through cells even when both segment endpoints are outside.
    for start, end in zip(vertices[:-1], vertices[1:]):
        abc = np.einsum("cij,cj->ci", transforms, start - cells[:, 0])
        delta = np.einsum("cij,j->ci", transforms, end - start)
        first = np.column_stack((1 - abc.sum(axis=1), abc))
        slope = np.column_stack((-delta.sum(axis=1), delta))
        low, high = np.zeros(len(cells)), np.ones(len(cells))
        possible = np.ones(len(cells), dtype=bool)
        for coordinate in range(4):
            value, change = first[:, coordinate], slope[:, coordinate]
            parallel = np.abs(change) < 1e-14
            possible &= ~(parallel & (value < -1e-12))
            crossing = np.divide(
                -1e-12 - value, change, out=np.zeros_like(value), where=~parallel
            )
            low = np.where(change > 1e-14, np.maximum(low, crossing), low)
            high = np.where(change < -1e-14, np.minimum(high, crossing), high)
        if np.any(possible & (low <= high)):
            raise ValueError("loop touches or intersects the conductor")
    integrated = np.empty((len(cells), 3))
    for first in range(0, len(cells), cell_chunk_size):
        last = min(first + cell_chunk_size, len(cells))
        local = cells[first:last]
        points = local[:, :1] + np.einsum(
            "qv,cvi->cqi", barycentric, local - local[:, :1]
        )
        potential = vector_potential(vertices, points.reshape(-1, 3)).reshape(
            points.shape
        )
        # The existing line kernel uses the historical mu0/(4*pi)=1e-7.
        potential *= mu_0 / (4 * np.pi * 1e-7)
        integrated[first:last] = np.einsum(
            "q,cqi,c->ci", weights, potential, mesh.volumes[first:last]
        )
    return sum(
        operator.T @ integrated[:, k] for k, operator in enumerate(basis.operators)
    )


def modal_response(
    rates: ArrayLike,
    modal_coupling: ArrayLike,
    current_derivative: Callable[[float], float],
    times: ArrayLike,
    *,
    initial: ArrayLike | None = None,
    rtol: float = 1e-8,
    atol: float = 1e-12,
    max_step: float = np.inf,
) -> NDArray[np.float64]:
    """Integrate da/dt + rates*a = -modal_coupling*dI/dt.

    modal_coupling = V.T @ m for M-orthonormal generalized eigenvectors V.
    Output shape is (nmodes, ntimes); coefficient currents are V @ output.
    The initial amplitudes apply at times[0], and default to zero. Specify
    max_step to resolve short pulses or sharp waveform features.
    """
    for value in (rates, modal_coupling, times, initial):
        if value is not None and np.iscomplexobj(value):
            raise ValueError("modal inputs must be real")
    rates = np.asarray(rates, dtype=float)
    coupling = np.asarray(modal_coupling, dtype=float)
    times = np.asarray(times, dtype=float)
    if (
        rates.ndim != 1
        or not rates.size
        or not np.isfinite(rates).all()
        or np.any(rates <= 0)
    ):
        raise ValueError("rates must be a nonempty positive finite vector")
    if coupling.shape != rates.shape or not np.isfinite(coupling).all():
        raise ValueError("modal_coupling must be finite and match rates")
    if (
        times.ndim != 1
        or len(times) < 2
        or not np.isfinite(times).all()
        or np.any(np.diff(times) <= 0)
    ):
        raise ValueError("times must be finite and strictly increasing")
    initial = (
        np.zeros_like(rates) if initial is None else np.asarray(initial, dtype=float)
    )
    if initial.shape != rates.shape or not np.isfinite(initial).all():
        raise ValueError("initial must be finite and match rates")
    solution = solve_ivp(
        lambda t, a: -rates * a - coupling * current_derivative(t),
        (times[0], times[-1]),
        initial,
        t_eval=times,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y
