"""Geometry-only Coulomb interactions of uniformly filled tetrahedra.

Near-source point potentials use analytic triangle integration. Outer volume
integration uses positive Gauss–Jacobi product quadrature with an empirical
order-refinement check. No coincident 1/r samples or diagonal replacements are
used. The error estimate is a convergence diagnostic, not a rigorous bound.
"""

from __future__ import annotations

from collections.abc import Iterator
from functools import lru_cache

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import roots_jacobi

from ..integrals_cached import TriangleIntegrals
from .mesh import TETRA_FACES, TetraMesh


def _positive_integer(value: int, name: str) -> None:
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        or value < 1
    ):
        raise ValueError(f"{name} must be a positive integer")


@lru_cache(maxsize=16)
def tetrahedron_quadrature(
    order: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Positive interior barycentric points and weights summing to one.

    There are order**3 points. The collapsed coordinates are
    (1-a, a(1-b), ab(1-c), abc), with Jacobian a²b. Gauss–Jacobi rules absorb
    those factors, so polynomials of total degree <= 2*order-1 are integrated
    exactly. The rule depends on vertex ordering; converged integrals do not.
    Returned cached arrays are read-only.
    """
    _positive_integer(order, "order")
    nodes, weights = [], []
    for beta in (2, 1, 0):
        x, w = roots_jacobi(order, 0, beta)
        nodes.append((x + 1) / 2)
        weights.append(w / 2 ** (beta + 1))
    a, b, c = np.meshgrid(*nodes, indexing="ij")
    barycentric = np.stack(
        (1 - a, a * (1 - b), a * b * (1 - c), a * b * c), axis=-1
    ).reshape(-1, 4)
    weights = 6 * np.einsum("i,j,k->ijk", *weights).ravel()
    barycentric.flags.writeable = False
    weights.flags.writeable = False
    return barycentric, weights


def _cell_indices(mesh: TetraMesh, cells: ArrayLike | None) -> NDArray[np.int64]:
    if cells is None:
        return np.arange(len(mesh.tetrahedra), dtype=np.int64)
    cells = np.asarray(cells)
    if cells.ndim != 1 or (cells.size and not np.issubdtype(cells.dtype, np.integer)):
        raise ValueError("cell indices must be a one-dimensional integer array")
    if np.any(cells < 0) or np.any(cells >= len(mesh.tetrahedra)):
        raise ValueError("cell index is out of range")
    return cells.astype(np.int64, copy=False)


def tetrahedron_potential(
    mesh: TetraMesh,
    points: ArrayLike,
    source_cells: ArrayLike | None = None,
    *,
    point_chunk_size: int = 256,
    source_chunk_size: int = 16,
) -> NDArray[np.float64]:
    """Uniform-source potentials integral_U 1/|x-y| dV_y, in squared length units.

    points has shape (npoints, 3); output has shape (npoints, nsources), following
    source_cells order (all cells by default). Interior and boundary points are
    allowed. Outward face heights give u_U(x) = 1/2 sum_f d_f(x) S_f(x).

    Each source is translated and scaled before calling the triangle routine to
    control its fixed logarithm regularization. Beyond eight source radii, a
    regular 64-point source quadrature avoids cancellation in the signed face
    sum. Its smooth-kernel error is much smaller than the default outer tolerance.
    Points and sources are chunked to bound displacement-array memory.
    """
    _positive_integer(point_chunk_size, "point_chunk_size")
    _positive_integer(source_chunk_size, "source_chunk_size")
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("points must be a finite array of shape (npoints, 3)")
    sources = _cell_indices(mesh, source_cells)
    result = np.empty((len(points), len(sources)))
    far_barycentric, far_weights = tetrahedron_quadrature(4)
    vertices = mesh.vertices[mesh.tetrahedra[sources]]
    centers = vertices.mean(axis=1)
    scales = np.linalg.norm(vertices - centers[:, None], axis=-1).max(axis=1)
    vertices = (vertices - centers[:, None]) / scales[:, None, None]
    # Advanced indexing otherwise puts the face axes before the source axis in
    # memory. A contiguous layout lets displacement reshapes remain views.
    faces = np.ascontiguousarray(vertices[:, np.asarray(TETRA_FACES)])
    # Opposite-edge order matches roll(R, 2) - roll(R, 1) for R = x - vertex.
    # Source geometry is shared by every point chunk, including the normals.
    edges = np.empty_like(faces)
    edges[:, :, 0] = faces[:, :, 2] - faces[:, :, 1]
    edges[:, :, 1] = faces[:, :, 0] - faces[:, :, 2]
    edges[:, :, 2] = faces[:, :, 1] - faces[:, :, 0]
    normals = np.cross(edges[:, :, 2], -edges[:, :, 1])
    normals /= np.linalg.norm(normals, axis=-1)[..., None]
    far_points = np.einsum("qv,svi->sqi", far_barycentric, vertices)
    volumes = mesh.volumes[sources]
    for column in range(0, len(sources), source_chunk_size):
        columns = slice(column, min(column + source_chunk_size, len(sources)))
        for start in range(0, len(points), point_chunk_size):
            stop = min(start + point_chunk_size, len(points))
            x = (points[start:stop, None] - centers[None, columns]) / scales[
                None, columns, None
            ]
            far = np.linalg.norm(x, axis=-1) > 8
            values = np.empty(far.shape)
            p, s = np.nonzero(~far)
            if p.size:
                if not np.any(far):
                    # Keep the evaluation axis: triangle edges and normals are
                    # shared across points, avoiding repeated geometry work.
                    displacements = (
                        x[:, :, None, None, :] - faces[None, columns]
                    ).reshape(len(x), -1, 3, 3)
                    triangle_normals = normals[columns].reshape(-1, 3)
                    triangle_edges = edges[columns].reshape(-1, 3, 3)
                else:
                    # Pack selected point/source/face combinations into one
                    # triangle axis, matching normals and cached source edges.
                    displacements = (
                        x[p, s, None, None, :] - faces[columns][s]
                    ).reshape(1, -1, 3, 3)
                    triangle_normals = normals[columns][s].reshape(-1, 3)
                    triangle_edges = edges[columns][s].reshape(-1, 3, 3)
                triangle = TriangleIntegrals(
                    displacements, triangle_normals, edges=triangle_edges
                )
                heights = -triangle.d_distance.reshape(-1, 4)
                triangle_potentials = triangle.triangle_potential_uniform.reshape(-1, 4)
                # Zero heights suppress boundary solid-angle branch ambiguities.
                values[p, s] = (
                    0.5
                    * np.sum(heights * triangle_potentials, axis=-1)
                    * scales[columns][s] ** 2
                )
                # Retain caches only for this evaluation batch.
                del triangle
            p, s = np.nonzero(far)
            if p.size:
                distances = np.linalg.norm(
                    x[p, s, None, :] - far_points[columns][s], axis=-1
                )
                values[p, s] = (
                    volumes[columns][s]
                    / scales[columns][s]
                    * ((1 / distances) @ far_weights)
                )
            if not np.all(np.isfinite(values)) or np.any(values <= 0):
                raise FloatingPointError(
                    "nonpositive or nonfinite tetrahedron potential"
                )
            result[start:stop, columns] = values
    return result


def _integrate_sources(
    mesh: TetraMesh, rows: NDArray, columns: NDArray, order: int
) -> NDArray[np.float64]:
    barycentric, weights = tetrahedron_quadrature(order)
    # Work relative to a cell vertex so a large coordinate offset does not enter
    # the barycentric weighted sum multiple times.
    vertices = mesh.vertices[mesh.tetrahedra[rows]]
    points = vertices[:, :1] + np.einsum(
        "qv,tvi->tqi", barycentric, vertices - vertices[:, :1]
    )
    potentials = tetrahedron_potential(mesh, points.reshape(-1, 3), columns)
    potentials = potentials.reshape(len(rows), len(weights), len(columns))
    return mesh.volumes[rows, None] * np.einsum("q,tqs->ts", weights, potentials)


def _directed_interactions(
    mesh: TetraMesh,
    rows: NDArray,
    columns: NDArray,
    order: int,
    rtol: float,
    max_order: int,
) -> NDArray[np.float64]:
    previous = _integrate_sources(mesh, rows, columns, order)
    current = _integrate_sources(mesh, rows, columns, order + 2)
    relative_change = np.abs(current - previous) / np.abs(current)
    active = relative_change > rtol
    for refined_order in range(order + 4, max_order + 1, 2):
        if not np.any(active):
            break
        for row in np.flatnonzero(np.any(active, axis=1)):
            selected = np.flatnonzero(active[row])
            values = _integrate_sources(
                mesh, rows[row : row + 1], columns[selected], refined_order
            )[0]
            relative_change[row, selected] = np.abs(
                values - current[row, selected]
            ) / np.abs(values)
            current[row, selected] = values
        active = relative_change > rtol
    if np.any(active):
        i, j = np.argwhere(active)[0]
        raise RuntimeError(
            f"tetrahedron interaction ({rows[i]}, {columns[j]}) did not converge "
            f"at max_order={max_order}: relative change {relative_change[i,j]:.3g} "
            f"exceeds rtol={rtol:g}; increase max_order or relax rtol"
        )
    return current


def tetrahedron_interaction_block(
    mesh: TetraMesh,
    rows: ArrayLike,
    columns: ArrayLike,
    *,
    quadrature_order: int = 4,
    rtol: float = 1e-4,
    max_order: int = 12,
) -> NDArray[np.float64]:
    """Symmetric-by-reciprocity cell interactions K_TU, in length**5 units.

    Analytic inner integration is followed by increasing outer quadrature order
    by two until successive results change by at most rtol, for each directed
    cell pair. Both directions are evaluated and averaged. Failure to converge
    raises RuntimeError; no clipping, diagonal replacement, or regularization of
    magnetic energy is applied. The refinement estimate is empirical.
    """
    _positive_integer(quadrature_order, "quadrature_order")
    _positive_integer(max_order, "max_order")
    if max_order < quadrature_order + 2:
        raise ValueError("max_order must be at least quadrature_order + 2")
    if not np.isfinite(rtol) or rtol <= 0:
        raise ValueError("rtol must be positive and finite")
    rows, columns = _cell_indices(mesh, rows), _cell_indices(mesh, columns)
    if not len(rows) or not len(columns):
        return np.empty((len(rows), len(columns)))
    forward = _directed_interactions(
        mesh, rows, columns, quadrature_order, rtol, max_order
    )
    if np.array_equal(rows, columns):
        reverse = forward.T
    else:
        reverse = _directed_interactions(
            mesh, columns, rows, quadrature_order, rtol, max_order
        ).T
    return 0.5 * (forward + reverse)


def iter_tetrahedron_interactions(
    mesh: TetraMesh,
    *,
    block_size: int = 16,
    quadrature_order: int = 4,
    rtol: float = 1e-4,
    max_order: int = 12,
) -> Iterator[tuple[slice, slice, NDArray[np.float64]]]:
    """Yield upper-triangular cell blocks (row slice, column slice, K block).

    Off-diagonal blocks represent both orientations via transposition; consumers
    must include both. No complete cell-pair matrix is retained.
    """
    _positive_integer(block_size, "block_size")
    ncells = len(mesh.tetrahedra)
    for start in range(0, ncells, block_size):
        rows = slice(start, min(start + block_size, ncells))
        for other in range(start, ncells, block_size):
            columns = slice(other, min(other + block_size, ncells))
            block = tetrahedron_interaction_block(
                mesh,
                np.arange(rows.start, rows.stop),
                np.arange(columns.start, columns.stop),
                quadrature_order=quadrature_order,
                rtol=rtol,
                max_order=max_order,
            )
            yield rows, columns, block


def tetrahedron_interaction_matrix(
    mesh: TetraMesh,
    *,
    block_size: int = 16,
    quadrature_order: int = 4,
    rtol: float = 1e-4,
    max_order: int = 12,
) -> NDArray[np.float64]:
    """Materialize the dense geometry-only K matrix; use the iterator to stream."""
    matrix = np.empty((len(mesh.tetrahedra), len(mesh.tetrahedra)))
    for rows, columns, block in iter_tetrahedron_interactions(
        mesh,
        block_size=block_size,
        quadrature_order=quadrature_order,
        rtol=rtol,
        max_order=max_order,
    ):
        matrix[rows, columns] = block
        if rows.start != columns.start:
            matrix[columns, rows] = block.T
    return matrix
