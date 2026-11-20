"""Cached analytic triangle integrals, derived from :mod:`bfieldtools.integrals`.

The functional module remains independent. Each TriangleIntegrals instance owns
one evaluation batch and shares its geometry, distances, and projections across
integral formulas. Release it after consuming the batch to release its caches.
API consolidation is deferred; see TRIANGLE_INTEGRATION_API.md in this package.
"""

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = ["TriangleIntegrals"]


def determinant(a):
    """Faster determinant for the two last dimensions of 'a'"""
    det = a[..., 0, 0] * (a[..., 1, 1] * a[..., 2, 2] - a[..., 2, 1] * a[..., 1, 2])
    det += a[..., 0, 1] * (a[..., 1, 2] * a[..., 2, 0] - a[..., 2, 2] * a[..., 1, 0])
    det += a[..., 0, 2] * (a[..., 1, 0] * a[..., 2, 1] - a[..., 2, 0] * a[..., 1, 1])
    return det


def norm(vecs):
    """Faster vector norm for the last dimension of 'vecs'"""
    return np.sqrt(np.einsum("...i,...i", vecs, vecs))


def cross(r1, r2):
    """Cross product without overhead for the last dimensions of 'r1' and 'r2'"""
    result = np.zeros(r1.shape)
    result[..., 0] = r1[..., 1] * r2[..., 2] - r1[..., 2] * r2[..., 1]
    result[..., 1] = r1[..., 2] * r2[..., 0] - r1[..., 0] * r2[..., 2]
    result[..., 2] = r1[..., 0] * r2[..., 1] - r1[..., 1] * r2[..., 0]
    return result


@dataclass(frozen=True, eq=False, init=False)
class TriangleIntegrals:
    """Lazy integrals for fixed displacements R and triangle normals tn.

    R has shape (Neval, Ntri, 3, 3), or (Neval, 3, 3) for a single triangle.
    tn has shape (Ntri, 3), or (3,) for a single triangle. As in the functional
    implementation, source geometry is shared along the first evaluation axis.

    Optional edges seed the cached ``edges`` property. Their shape is R.shape[1:]
    and orientation is roll(R[0], 2, -2) - roll(R[0], 1, -2). Optional triangle
    areas ta similarly seed the ``ta`` property; otherwise areas are derived
    from the edges. Cached ``x_distance`` is unnormalized, and
    ``normalized_x_distance`` includes division by twice the triangle area.

    Inputs are borrowed without copying or changing their writeability. They
    and cached results must remain unchanged for the lifetime of this instance.
    Frozen fields prevent accidental replacement after caches are populated.
    ``reg`` and ``symmetrize`` retain gamma0's functional defaults; ``planar``
    omits the solid-angle term in the uniform triangle potential.

    For the vertex-dipole approximation only, R may instead have shape
    (Neval, Nvertex, 3), with vertex normals tn and supplied vertex areas ta.
    """

    R: NDArray = field(repr=False)
    tn: NDArray = field(repr=False)
    reg: float = 1e-13
    symmetrize: bool = True
    planar: bool = False
    approx_reg: float = 1e-12

    def __init__(
        self,
        R: ArrayLike,
        tn: ArrayLike,
        *,
        edges: ArrayLike | None = None,
        ta: ArrayLike | None = None,
        reg: float = 1e-13,
        symmetrize: bool = True,
        planar: bool = False,
        approx_reg: float = 1e-12,
    ):
        object.__setattr__(self, "R", np.asarray(R))
        object.__setattr__(self, "tn", np.asarray(tn))
        object.__setattr__(self, "reg", reg)
        object.__setattr__(self, "symmetrize", symmetrize)
        object.__setattr__(self, "planar", planar)
        object.__setattr__(self, "approx_reg", approx_reg)
        # Seed cached_property storage without introducing a second public edge
        # representation or allocating copies of caller-provided geometry.
        if edges is not None:
            self.__dict__["edges"] = np.asarray(edges)
        if ta is not None:
            self.__dict__["ta"] = np.asarray(ta)

    @cached_property
    def edges(self):
        """Opposite-edge vectors, shared along the evaluation axis."""
        return np.roll(self.R[0], 2, -2) - np.roll(self.R[0], 1, -2)

    @cached_property
    def edge_lengths(self):
        return norm(self.edges)

    @cached_property
    def ta(self):
        """Triangle areas; supplied areas take precedence over this calculation."""
        return norm(cross(self.edges[..., 2, :], -self.edges[..., 1, :])) / 2

    @cached_property
    def vertex_distances(self):
        """Distances reused by both edge potentials and solid angles."""
        return norm(self.R)

    @cached_property
    def edge_dotprods1(self):
        return np.roll(
            np.einsum("...i,...i", self.R, np.roll(self.edges, -1, -2)), 1, -1
        )

    @cached_property
    def edge_dotprods2(self):
        return np.roll(
            np.einsum("...i,...i", self.R, np.roll(self.edges, -2, -2)), 2, -1
        )

    @cached_property
    def scaled_distances1(self):
        return np.roll(self.vertex_distances, 1, -1) * self.edge_lengths

    @cached_property
    def scaled_distances2(self):
        return np.roll(self.vertex_distances, 2, -1) * self.edge_lengths

    @cached_property
    def gamma0(self):
        """Unit-parameter edge potentials with the original regularization."""
        nn1, nn2 = self.scaled_distances1, self.scaled_distances2
        dotprods1, dotprods2 = self.edge_dotprods1, self.edge_dotprods2
        res = np.log((nn1 + dotprods1 + self.reg) / (nn2 + dotprods2 + self.reg))
        if self.symmetrize:
            mask = (np.abs(dotprods1 + nn1) < 1e-12) * (dotprods1 + dotprods2 < 0)
            res[mask] = -np.log(
                (nn1[mask] - dotprods1[mask] + self.reg)
                / (nn2[mask] - dotprods2[mask] + self.reg)
            )
        res /= self.edge_lengths
        return -res

    @cached_property
    def determinant(self):
        return determinant(self.R)

    @cached_property
    def omega(self):
        """Solid angles, sharing vertex distances with gamma0."""
        d = self.vertex_distances
        denom = np.prod(d, axis=-1)
        for i in range(3):
            j, k = (i + 1) % 3, (i + 2) % 3
            denom += np.einsum(
                "...i,...i,...", self.R[..., i, :], self.R[..., j, :], d[..., k]
            )
        return -2 * np.arctan2(self.determinant, denom)

    @cached_property
    def edge_conormals(self):
        # cross allocates its result, preserving supplied read-only edges.
        return -cross(self.edges, self.tn[..., None, :])

    @cached_property
    def x_distance(self):
        """Unnormalized projections onto the triangle edge conormals."""
        return np.roll(
            np.einsum("...i,...i", self.R, np.roll(self.edge_conormals, -1, -2)),
            1,
            -1,
        )

    @cached_property
    def normalized_x_distance(self):
        return self.x_distance / (2 * self.ta[..., None])

    @cached_property
    def d_distance(self):
        """Signed face-plane distances, also used for tetrahedron face heights."""
        return np.einsum("...i,...i->...", self.R[..., 0, :], self.tn)

    @cached_property
    def c_coeffs(self):
        return np.einsum(
            "...ik,...jk->...ij",
            self.edges,
            self.edges / (2 * self.ta[..., None, None]),
        )

    @cached_property
    def triangle_potential_uniform(self):
        result = np.einsum("...i,...i", self.gamma0, self.x_distance)
        if not self.planar:
            result += self.d_distance * self.omega
        return result

    @cached_property
    def triangle_centers(self):
        """Displacements to triangle centroids."""
        return self.R.mean(axis=-2)

    @cached_property
    def triangle_potential_approx(self):
        return self.ta / (norm(self.triangle_centers) + self.approx_reg)

    @cached_property
    def dipole_quadrature_displacements(self):
        weights = np.array([[0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5]])
        return np.einsum("...ij,ik->...kj", self.R, weights)

    @cached_property
    def potential_dipoles(self):
        Rquad = self.dipole_quadrature_displacements
        pot = np.einsum("...i,...ji->...j", self.tn, Rquad) / norm(Rquad) ** 3
        return pot * (self.ta[..., None] / 3)

    @cached_property
    def potential_vertex_dipoles(self):
        pot = np.einsum("ik,lik->li", self.tn, self.R) / self.vertex_distances**3
        return pot * self.ta

    @cached_property
    def triangle_potential_dipole_linear(self):
        result = np.einsum(
            "...i,...ij,...->...j",
            self.gamma0,
            self.c_coeffs,
            self.d_distance,
            optimize=True,
        )
        result -= self.normalized_x_distance * self.omega[..., :, None]
        return result
