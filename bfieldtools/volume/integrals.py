"""Fixed-order Coulomb interactions for tetrahedral volume currents.

Self and shared-vertex neighbors use analytic source potentials with fixed outer
quadrature. Other pairs use direct double volume quadrature. The defaults target
slow, smooth eddy-current modes; there is no adaptive accuracy estimator.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from functools import cached_property, lru_cache

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


@dataclass(frozen=True)
class _SourceGeometry:
    """Normalized source geometry shared by all potential evaluations.

    take() selects matching source arrays for a cell block or neighbor group.
    """

    centers: NDArray
    scales: NDArray
    faces: NDArray
    edges: NDArray
    normals: NDArray
    vertices: NDArray
    volumes: NDArray

    @classmethod
    def from_mesh(cls, mesh: TetraMesh, sources: NDArray) -> _SourceGeometry:
        vertices = mesh.vertices[mesh.tetrahedra[sources]]
        centers = mesh.cell_centers[sources]
        scales = mesh.cell_radii[sources]
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
        volumes = mesh.volumes[sources]
        return cls(centers, scales, faces, edges, normals, vertices, volumes)

    @cached_property
    def far_points(self):
        """The stable point-potential fallback is prepared only when needed."""
        barycentric, _ = tetrahedron_quadrature(4)
        return np.einsum("qv,svi->sqi", barycentric, self.vertices)

    def take(self, indices: NDArray) -> _SourceGeometry:
        return _SourceGeometry(
            self.centers[indices],
            self.scales[indices],
            self.faces[indices],
            self.edges[indices],
            self.normals[indices],
            self.vertices[indices],
            self.volumes[indices],
        )


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
    sum. This fixed fallback controls cancellation; it is independent of the
    cell-interaction quadrature orders.
    Points and sources are chunked to bound displacement-array memory.
    """
    _positive_integer(point_chunk_size, "point_chunk_size")
    _positive_integer(source_chunk_size, "source_chunk_size")
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("points must be a finite array of shape (npoints, 3)")
    sources = _cell_indices(mesh, source_cells)
    geometry = _SourceGeometry.from_mesh(mesh, sources)
    return _evaluate_potential(geometry, points, point_chunk_size, source_chunk_size)


def _evaluate_potential(
    geometry: _SourceGeometry,
    points: NDArray,
    point_chunk_size: int = 256,
    source_chunk_size: int = 16,
) -> NDArray[np.float64]:
    centers, scales = geometry.centers, geometry.scales
    faces, edges, normals = geometry.faces, geometry.edges, geometry.normals
    volumes = geometry.volumes
    nsources = len(centers)
    _, far_weights = tetrahedron_quadrature(4)
    result = np.empty((len(points), nsources))
    for column in range(0, nsources, source_chunk_size):
        columns = slice(column, min(column + source_chunk_size, nsources))
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
                    x[p, s, None, :] - geometry.far_points[columns][s], axis=-1
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


@dataclass(frozen=True)
class _CellIntegrals:
    """Prepared rules and source geometry, reused across streamed cell blocks."""

    tetrahedra: NDArray
    sources: _SourceGeometry
    regular_offsets: NDArray
    regular_weights: NDArray
    touching_points: NDArray
    touching_weights: NDArray

    @classmethod
    def from_mesh(cls, mesh, cells=None, *, quadrature_order=2, touching_order=4):
        _positive_integer(quadrature_order, "quadrature_order")
        _positive_integer(touching_order, "touching_order")
        cells = _cell_indices(mesh, cells)
        vertices = mesh.vertices[mesh.tetrahedra[cells]]
        sources = _SourceGeometry.from_mesh(mesh, cells)

        def rule(order):
            barycentric, weights = tetrahedron_quadrature(order)
            offsets = np.einsum(
                "qv,tvi->tqi", barycentric, vertices - sources.centers[:, None]
            )
            return offsets, weights

        regular, regular_weights = rule(quadrature_order)
        touching, touching_weights = rule(touching_order)
        return cls(
            mesh.tetrahedra[cells],
            sources,
            regular,
            regular_weights,
            sources.centers[:, None] + touching,
            touching_weights,
        )

    def _regular(self, rows, columns):
        """Double quadrature, with bounded temporary point-pair arrays."""
        result = np.empty((len(rows), len(columns)))
        centers, volumes = self.sources.centers, self.sources.volumes
        for i, row in enumerate(rows):
            for start in range(0, len(columns), 16):
                batch = slice(start, start + 16)
                selected = columns[batch]
                displacement = (
                    (centers[row] - centers[selected])[:, None, None, :]
                    + self.regular_offsets[row, None, :, None, :]
                    - self.regular_offsets[selected, None, :, :]
                )
                result[i, batch] = (
                    volumes[row]
                    * volumes[selected]
                    * np.einsum(
                        "q,sqp,p->s",
                        self.regular_weights,
                        1 / np.linalg.norm(displacement, axis=-1),
                        self.regular_weights,
                    )
                )
        return result

    def _analytic(self, rows, columns):
        """Fixed outer quadrature of analytic source potentials; average directions."""

        def directed(testing, sources):
            values = _evaluate_potential(
                self.sources.take(sources), self.touching_points[testing].reshape(-1, 3)
            ).reshape(len(testing), len(self.touching_weights), len(sources))
            return self.sources.volumes[testing, None] * np.einsum(
                "q,tqs->ts", self.touching_weights, values
            )

        forward = directed(rows, columns)
        reverse = (
            forward.T if np.array_equal(rows, columns) else directed(columns, rows).T
        )
        return 0.5 * (forward + reverse)

    def block(self, rows, columns):
        if not len(rows) or not len(columns):
            return np.empty((len(rows), len(columns)))
        touching = np.any(
            self.tetrahedra[rows, None, :, None]
            == self.tetrahedra[None, columns, None, :],
            axis=(2, 3),
        )
        if not np.any(touching):
            result = self._regular(rows, columns)
        elif np.all(touching):
            result = self._analytic(rows, columns)
        else:
            result = np.empty(touching.shape)
            for i in range(len(rows)):
                for mask, integrate in (
                    (touching[i], self._analytic),
                    (~touching[i], self._regular),
                ):
                    selected = columns[mask]
                    if len(selected):
                        result[i, mask] = integrate(rows[i : i + 1], selected)[0]
        if not np.all(np.isfinite(result)) or np.any(result <= 0):
            raise FloatingPointError("nonpositive or nonfinite tetrahedron interaction")
        return 0.5 * (result + result.T) if np.array_equal(rows, columns) else result


def tetrahedron_interaction_block(
    mesh: TetraMesh,
    rows: ArrayLike,
    columns: ArrayLike,
    *,
    quadrature_order: int = 2,
    touching_order: int = 4,
) -> NDArray[np.float64]:
    """Cell interactions K_TU in length**5 units, symmetric by reciprocity.

    Self and cells sharing any vertex use analytic source potentials integrated
    with a fixed order-touching_order outer rule. All other pairs use fixed
    double volume quadrature at quadrature_order (order**3 points per cell).
    Defaults are two and four, intended for slow, smooth current patterns.
    Increase either order explicitly to check sensitivity; no error bound,
    runtime refinement, diagonal replacement, or energy correction is applied.
    Supplied meshes must be geometrically conforming with shared vertex indices.
    """
    rows, columns = _cell_indices(mesh, rows), _cell_indices(mesh, columns)
    cells = np.unique(np.concatenate((rows, columns)))
    geometry = _CellIntegrals.from_mesh(
        mesh, cells, quadrature_order=quadrature_order, touching_order=touching_order
    )
    return geometry.block(np.searchsorted(cells, rows), np.searchsorted(cells, columns))


def iter_tetrahedron_interactions(
    mesh: TetraMesh,
    *,
    block_size: int = 16,
    quadrature_order: int = 2,
    touching_order: int = 4,
    spatial_order: bool = False,
) -> Iterator[tuple[slice | NDArray, slice | NDArray, NDArray[np.float64]]]:
    """Yield upper-triangular cell blocks (rows, columns, K block).

    Off-diagonal blocks represent both orientations via transposition. Source
    geometry and both fixed rules are prepared once; no full cell-pair matrix
    is retained. The default selectors are consecutive slices. spatial_order=True
    groups nearby centroids and yields original-numbering index arrays instead;
    scatter those blocks with np.ix_(rows, columns).
    """
    _positive_integer(block_size, "block_size")
    geometry = _CellIntegrals.from_mesh(
        mesh, quadrature_order=quadrature_order, touching_order=touching_order
    )
    ncells = len(mesh.tetrahedra)
    order = mesh.spatial_cell_order if spatial_order else np.arange(ncells)
    for start in range(0, ncells, block_size):
        rows = slice(start, min(start + block_size, ncells))
        for other in range(start, ncells, block_size):
            columns = slice(other, min(other + block_size, ncells))
            block = geometry.block(order[rows], order[columns])
            if spatial_order:
                yield order[rows], order[columns], block
            else:
                yield rows, columns, block


def tetrahedron_interaction_matrix(
    mesh: TetraMesh,
    *,
    block_size: int = 16,
    quadrature_order: int = 2,
    touching_order: int = 4,
    spatial_order: bool = True,
) -> NDArray[np.float64]:
    """Dense K in original cell order; use the iterator to stream.

    See tetrahedron_interaction_block() for the two fixed quadrature orders.
    """
    matrix = np.empty((len(mesh.tetrahedra), len(mesh.tetrahedra)))
    for rows, columns, block in iter_tetrahedron_interactions(
        mesh,
        block_size=block_size,
        quadrature_order=quadrature_order,
        touching_order=touching_order,
        spatial_order=spatial_order,
    ):
        selectors = np.ix_(rows, columns) if spatial_order else (rows, columns)
        matrix[selectors] = block
        if not np.array_equal(rows, columns):
            selectors = np.ix_(columns, rows) if spatial_order else (columns, rows)
            matrix[selectors] = block.T
    return matrix
