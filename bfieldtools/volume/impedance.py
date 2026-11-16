"""Material matrices for conserved, piecewise-constant volume currents."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.constants import mu_0
from scipy.sparse import csr_matrix

from .basis import CellCurrentBasis
from .integrals import iter_tetrahedron_interactions


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


def inductance_matrix(
    basis: CellCurrentBasis,
    *,
    interactions: ArrayLike | None = None,
    block_size: int = 16,
    quadrature_order: int = 4,
    rtol: float = 1e-4,
    max_order: int = 12,
) -> NDArray[np.float64]:
    """Dense free-space inductance M = mu_0/(4*pi) sum_k Bk.T K Bk, in H.

    By default, geometry-only cell interactions are integrated in blocks without
    storing the complete K matrix. quadrature_order, rtol, and max_order control
    the outer integration; see tetrahedron_interaction_block(). block_size bounds
    the cell blocks and dense contraction intermediates. The final basis matrix
    is dense and occupies 8*n_dofs**2 bytes.

    interactions optionally supplies a precomputed real symmetric (ncells,ncells)
    K matrix in the basis mesh's exact cell ordering. K includes both volume
    integrations and has units m**5. Numerical integration settings are unused
    when K is supplied. Positive magnetic energy still depends on its accuracy;
    this function never alters eigenvalues to enforce definiteness.

    Geometry must be in metres, and the supplied basis must represent conserved
    currents for the eddy-current formulation. Energy is x.T @ M @ x / 2 for
    coefficients x in A. Material permeability is mu_0 (nonmagnetic free space).
    """
    if (
        isinstance(block_size, (bool, np.bool_))
        or not isinstance(block_size, (int, np.integer))
        or block_size < 1
    ):
        raise ValueError("block_size must be a positive integer")
    ncells = len(basis.mesh.tetrahedra)
    if interactions is None:
        blocks = iter_tetrahedron_interactions(
            basis.mesh,
            block_size=block_size,
            quadrature_order=quadrature_order,
            rtol=rtol,
            max_order=max_order,
        )
    else:
        interactions = np.asarray(interactions)
        if np.iscomplexobj(interactions) or interactions.shape != (ncells, ncells):
            raise ValueError("interactions must be a real (ncells, ncells) array")
        interactions = np.asarray(interactions, dtype=float)
        if not np.all(np.isfinite(interactions)) or not np.allclose(
            interactions, interactions.T, rtol=1e-12, atol=0
        ):
            raise ValueError("interactions must be finite and symmetric")
        blocks = (
            (
                slice(i, min(i + block_size, ncells)),
                slice(j, min(j + block_size, ncells)),
                interactions[i : i + block_size, j : j + block_size],
            )
            for i in range(0, ncells, block_size)
            for j in range(i, ncells, block_size)
        )
    matrix = np.zeros((basis.n_dofs, basis.n_dofs))
    if not basis.n_dofs:
        return matrix
    for rows, columns, block in blocks:
        for operator in basis.operators:
            left, right = operator[rows], operator[columns]
            active_rows = np.unique(left.indices)
            active_columns = np.unique(right.indices)
            if not len(active_rows) or not len(active_columns):
                continue
            # Restrict the dense contraction to coefficients supported in this
            # cell block. Edge bases typically activate only a few per cell;
            # avoid a full n_dofs² temporary for every geometric interaction.
            left = left[:, active_rows].toarray()
            right = right[:, active_columns].toarray()
            contribution = left.T @ block @ right
            matrix[np.ix_(active_rows, active_columns)] += contribution
            if rows.start != columns.start:
                matrix[np.ix_(active_columns, active_rows)] += contribution.T
    matrix *= mu_0 / (4 * np.pi)
    return (matrix + matrix.T) / 2
