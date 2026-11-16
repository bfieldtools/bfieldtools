"""Material matrices for conserved, piecewise-constant volume currents."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from scipy.sparse import csr_matrix

from .basis import CellCurrentBasis


def resistance_matrix(basis: CellCurrentBasis, resistivity: ArrayLike) -> csr_matrix:
    r"""Assemble exact sparse resistance for piecewise-constant currents.

    Parameters
    ----------
    basis : CellCurrentBasis
        Piecewise-constant current representation. Assembly is independent of
        its topology, gauge, or construction method.
    resistivity : float or array-like, shape (ncells,)
        Positive, finite scalar resistivity in ohm metres. A scalar applies to
        all cells; an array assigns one scalar per cell. Tensor resistivity is
        not supported yet.

    Returns
    -------
    scipy.sparse.csr_matrix
        Matrix R in ohms, with shape (n_dofs, n_dofs). Joule loss for real
        coefficients x in amperes is x.T @ R @ x, in watts.

    Notes
    -----
    With constant cell curls, R_ab = sum_T rho_T V_T j_a(T) dot j_b(T).
    Equivalently R = sum_k Bk.T diag(rho * V) Bk, using the basis's cached
    component operators. Assembly uses weighted Gram matrices and no quadrature.
    Ordinary R eigenpairs depend on coefficient normalization and tree gauge;
    eddy-current decay rates require the inductance matrix as the eigenproblem
    metric, R v = lambda M v.
    """
    mesh = basis.mesh
    resistivity = np.asarray(resistivity)
    if np.iscomplexobj(resistivity):
        raise ValueError("resistivity must be real, finite, and positive")
    resistivity = np.asarray(resistivity, dtype=float)
    if resistivity.ndim > 1 or (
        resistivity.ndim == 1 and resistivity.shape != mesh.volumes.shape
    ):
        raise ValueError("resistivity must be a scalar or have shape (ncells,)")
    if not np.all(np.isfinite(resistivity)) or np.any(resistivity <= 0):
        raise ValueError("resistivity must be real, finite, and positive")
    weights = np.sqrt(resistivity * mesh.volumes)
    matrix = csr_matrix((basis.n_dofs, basis.n_dofs), dtype=float)
    for operator in basis.operators:
        weighted = operator.multiply(weights[:, None]).tocsr()
        matrix += weighted.T @ weighted
    matrix.eliminate_zeros()
    return matrix.tocsr()
