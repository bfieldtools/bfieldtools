"""Analytic constant-face Laplace operators with Green function 1/(4*pi*r)."""

import numpy as np

from ..integrals import gamma0, omega, triangle_potential_uniform


def _geometry(mesh, points, Nchunks):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("points must be a finite array of shape (Npoints, 3).")
    triangles = np.asarray(mesh.triangles)
    if len(triangles) == 0 or not np.all(np.isfinite(triangles)):
        raise ValueError("mesh must contain finite, nondegenerate triangles.")
    if np.any(mesh.area_faces <= 0):
        raise ValueError("mesh must contain finite, nondegenerate triangles.")
    if Nchunks is None:
        # Limit displacement temporaries to about 100,000 target/source pairs.
        Nchunks = max(1, int(np.ceil(len(points) * len(triangles) / 100000)))
    if not isinstance(Nchunks, (int, np.integer)) or Nchunks < 1:
        raise ValueError("Nchunks must be a positive integer.")
    chunks = np.array_split(np.arange(len(points)), min(Nchunks, max(1, len(points))))
    # Normalize each source triangle: existing kernels use small absolute
    # regularizers, so evaluating at unit triangle scale preserves SI scaling.
    lengths = np.linalg.norm(triangles[:, 1] - triangles[:, 0], axis=-1)
    return points, triangles, lengths, chunks


def _gradient(r, normals):
    edges = np.roll(r[0], 2, axis=-2) - np.roll(r[0], 1, axis=-2)
    # Edges run counterclockwise; edge x n is the outward in-plane conormal.
    conormals = np.cross(edges, normals[:, None])
    return (
        omega(r)[..., None] * normals - np.einsum("pfe,fea->pfa", gamma0(r), conormals)
    ) / (4 * np.pi)


def _on_surface(r, normals):
    """Whether any target is on any source triangle, including edges."""
    a = r[:, :, 0]
    e1 = r[:, :, 0] - r[:, :, 1]
    e2 = r[:, :, 0] - r[:, :, 2]
    aa = np.einsum("pfi,pfi->pf", e1, e1)
    ab = np.einsum("pfi,pfi->pf", e1, e2)
    bb = np.einsum("pfi,pfi->pf", e2, e2)
    ap = np.einsum("pfi,pfi->pf", a, e1)
    bp = np.einsum("pfi,pfi->pf", a, e2)
    denom = aa * bb - ab**2
    u = (bb * ap - ab * bp) / denom
    v = (aa * bp - ab * ap) / denom
    tol = 64 * np.finfo(float).eps
    plane = abs(np.einsum("pfi,fi->pf", a, normals)) <= tol
    return np.any(plane & (u >= -tol) & (v >= -tol) & (u + v <= 1 + tol))


def single_layer_potential(mesh, points, Nchunks=None):
    """Map constant face density q [A/m] to potential u [A].

    Returns an array of shape (Npoints, Nfaces), integrating 1/(4*pi*r)
    over every source triangle. Source areas are included. Evaluation on a
    face interior is allowed. The mesh need not be closed.
    """
    points, triangles, lengths, chunks = _geometry(mesh, points, Nchunks)
    result = np.empty((len(points), len(triangles)))
    for ids in chunks:
        if len(ids):
            r = (points[ids, None, None] - triangles[None]) / lengths[
                None, :, None, None
            ]
            result[ids] = triangle_potential_uniform(r, mesh.face_normals) * (
                lengths / (4 * np.pi)
            )
    return result


def single_layer_gradient(mesh, points, Nchunks=None):
    """Map constant face density to grad(u) at off-surface points.

    Returns (Npoints, 3, Nfaces). The induced H field is the NEGATIVE
    contraction with q, since H = H_0 - grad(u). Reuses analytic edge
    integrals and solid angles, including for near-surface evaluations.
    Targets on a source triangle are rejected: one-sided traces require
    explicit jump terms and are not supplied by this function.
    """
    points, triangles, lengths, chunks = _geometry(mesh, points, Nchunks)
    result = np.empty((len(points), 3, len(triangles)))
    for ids in chunks:
        if len(ids):
            r = (points[ids, None, None] - triangles[None]) / lengths[
                None, :, None, None
            ]
            if _on_surface(r, mesh.face_normals):
                raise ValueError("Gradient targets must be off the source surface.")
            result[ids] = _gradient(r, mesh.face_normals).swapaxes(1, 2)
    return result


def adjoint_double_layer(boundary, Nchunks=None):
    """Assemble K' on a MagneticBoundary by face-centroid collocation.

    Returns (Nfaces, Nfaces), using metal-outward TARGET normals and the
    Green function 1/(4*pi*r). Its planar self terms are explicitly zero;
    the one-sided derivative jump +/-I/2 is not included. For finite
    permeability solve (I + 2*diag(beta)*K') q = 2*beta*g,
    where g = n dot H_0. No material factors are included in this operator.
    """
    mesh = boundary.mesh
    points, triangles, lengths, chunks = _geometry(mesh, mesh.triangles_center, Nchunks)
    result = np.empty((len(points), len(triangles)))
    for ids in chunks:
        if len(ids):
            r = (points[ids, None, None] - triangles[None]) / lengths[
                None, :, None, None
            ]
            result[ids] = np.einsum(
                "pfa,pa->pf", _gradient(r, mesh.face_normals), mesh.face_normals[ids]
            )
    np.fill_diagonal(result, 0.0)
    return result
