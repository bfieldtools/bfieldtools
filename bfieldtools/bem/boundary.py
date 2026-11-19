"""Validated, metal-outward boundary geometry for bodies in vacuum."""

import numpy as np
import trimesh

from ..integrals import omega


def _intersecting_triangles(triangle, others, tolerance):
    """Triangle separating-axis test, including coplanar pairs."""
    # Work near the origin to avoid cancellation in projections.
    others = others - triangle[0]
    triangle = triangle - triangle[0]
    edges = np.roll(triangle, -1, axis=0) - triangle
    other_edges = np.roll(others, -1, axis=1) - others
    normal = np.cross(edges[0], edges[1])
    other_normals = np.cross(other_edges[:, 0], other_edges[:, 1])
    count = len(others)
    axes = np.concatenate(
        [
            np.broadcast_to(normal, (count, 1, 3)),
            other_normals[:, None],
            np.cross(edges[None, :, None], other_edges[:, None]).reshape(count, 9, 3),
            np.broadcast_to(np.cross(normal, edges), (count, 3, 3)),
            np.cross(other_normals[:, None], other_edges),
        ],
        axis=1,
    )
    lengths = np.linalg.norm(axes, axis=-1)
    axes = axes / np.where(lengths > 0, lengths, 1)[..., None]
    a = np.einsum("vad,td->vat", axes, triangle)
    b = np.einsum("vad,vtd->vat", axes, others)
    separated = (a.max(axis=-1) < b.min(axis=-1) - tolerance) | (
        b.max(axis=-1) < a.min(axis=-1) - tolerance
    )
    return ~np.any(separated, axis=1)


def _validate_intersections(mesh, tolerance):
    """Reject crossings and contacts except ordinary mesh adjacency.

    Bounding boxes prune the quadratic search. Adjacent triangles are slightly
    contracted to exclude their permitted shared edge or vertex from the test.
    """
    triangles = mesh.triangles
    low, high = triangles.min(axis=1), triangles.max(axis=1)
    for i, triangle in enumerate(triangles[:-1]):
        candidates = (
            np.flatnonzero(
                np.all(high[i] >= low[i + 1 :] - tolerance, axis=1)
                & np.all(low[i] <= high[i + 1 :] + tolerance, axis=1)
            )
            + i
            + 1
        )
        if not len(candidates):
            continue
        shared = np.any(
            mesh.faces[i, None, :, None] == mesh.faces[candidates, None, :],
            axis=(1, 2),
        )
        for adjacent in (False, True):
            ids = candidates[shared == adjacent]
            if not len(ids):
                continue
            a, b = triangle, triangles[ids]
            if adjacent:
                a = a.mean(axis=0) + (a - a.mean(axis=0)) * (1 - 1e-8)
                centers = b.mean(axis=1, keepdims=True)
                b = centers + (b - centers) * (1 - 1e-8)
            hits = _intersecting_triangles(a, b, tolerance)
            if np.any(hits):
                raise ValueError(
                    "Boundary triangles intersect or touch: faces %d and %d."
                    % (i, ids[np.flatnonzero(hits)[0]])
                )


class MagneticBoundary:
    """A complete boundary of one or more disjoint metal bodies in vacuum.

    Parameters
    ----------
    meshes : trimesh.Trimesh or sequence of trimesh.Trimesh
        One mesh per connected metal body. A hollow body's mesh contains its
        outer surface and all cavity surfaces. Supply normals out of the metal,
        including inward-facing cavity surfaces. Faces are never reoriented.
        Surfaces must be watertight, nonintersecting and consistently wound.
    mu_r : float or array of shape (Nbodies,), optional
        Positive finite relative permeability, default 1000.

    Notes
    -----
    Input meshes are copied and concatenated without processing, preserving
    face order. ``body_ids`` maps each face to its input body; ``face_mu_r``
    and ``beta`` provide material values per face. Geometry validation uses
    floating-point tolerances relative to the boundary size. Treat the stored
    mesh as immutable; construct a new boundary after changing geometry.
    A cavity may contain another disjoint metal body, but overlapping metal
    volumes and shared interfaces are unsupported.
    """

    def __init__(self, meshes, mu_r=1000.0):
        if isinstance(meshes, trimesh.Trimesh):
            meshes = [meshes]
        else:
            meshes = list(meshes)
        if not meshes or not all(isinstance(m, trimesh.Trimesh) for m in meshes):
            raise ValueError("Provide a mesh or a nonempty sequence of meshes.")
        values = np.asarray(mu_r, dtype=float)
        if values.ndim == 0:
            values = np.full(len(meshes), values)
        if values.shape != (len(meshes),) or not np.all(
            np.isfinite(values) & (values > 0)
        ):
            raise ValueError(
                "mu_r must be positive and finite, with one value per body."
            )

        copies = []
        for mesh in meshes:
            if (
                len(mesh.faces) == 0
                or not np.all(np.isfinite(mesh.vertices))
                or np.any(mesh.faces < 0)
                or np.any(mesh.faces >= len(mesh.vertices))
            ):
                raise ValueError(
                    "Boundary meshes must have finite vertices and valid faces."
                )
            copy = trimesh.Trimesh(
                mesh.vertices.copy(), mesh.faces.copy(), process=False
            )
            if not copy.is_watertight or not copy.is_winding_consistent:
                raise ValueError(
                    "Boundary meshes must be watertight and consistently wound."
                )
            copies.append(copy)
        self.mesh = trimesh.util.concatenate(copies)
        self.body_ids = np.repeat(
            np.arange(len(copies)), [len(m.faces) for m in copies]
        )
        self.mu_r = values.copy()
        self.face_mu_r = self.mu_r[self.body_ids]
        # This form remains finite even near the largest representable mu_r.
        self.beta = 1 - 2 / (self.face_mu_r + 1)

        scale = np.max(np.ptp(self.mesh.vertices, axis=0))
        tolerance = 64 * np.finfo(float).eps * scale
        if scale == 0 or np.any(self.mesh.area_faces <= tolerance**2):
            raise ValueError("Boundary meshes must not contain degenerate triangles.")
        _validate_intersections(self.mesh, tolerance)

        components = trimesh.graph.connected_components(
            self.mesh.face_adjacency, nodes=np.arange(len(self.mesh.faces))
        )
        components = sorted(components, key=lambda c: int(np.min(c)))
        # After intersection checks, one representative vertex determines
        # containment of an entire connected surface. Solid angles avoid an
        # optional spatial-index dependency and do not depend on winding sign.
        inside = np.zeros((len(components), len(components)), dtype=bool)
        signs = []
        for i, faces in enumerate(components):
            triangles = self.mesh.triangles[faces]
            shifted = triangles - triangles[0, 0]
            volume = (
                np.einsum(
                    "ij,ij->i", shifted[:, 0], np.cross(shifted[:, 1], shifted[:, 2])
                ).sum()
                / 6
            )
            if abs(volume) <= tolerance**3:
                raise ValueError(
                    "Each closed boundary component must enclose a volume."
                )
            signs.append(np.sign(volume))
            for j, other in enumerate(components):
                if i != j:
                    r = triangles[0, 0] - self.mesh.triangles[other][None]
                    inside[i, j] = abs(omega(r).sum()) > 2 * np.pi
        depths = inside.sum(axis=1)
        if np.any(np.asarray(signs) != np.where(depths % 2 == 0, 1, -1)):
            raise ValueError(
                "Normals must point out of metal: outer surfaces outward, cavities inward; "
                "nested or overlapping metal volumes are unsupported."
            )
        component_bodies = np.array([self.body_ids[c[0]] for c in components])
        for body in range(len(copies)):
            if np.sum((component_bodies == body) & (np.asarray(signs) > 0)) != 1:
                raise ValueError("Supply one connected metal body per input mesh.")
        for i in np.flatnonzero(np.asarray(signs) < 0):
            parents = np.flatnonzero(inside[i] & (depths == depths[i] - 1))
            if len(parents) != 1 or component_bodies[parents[0]] != component_bodies[i]:
                raise ValueError(
                    "A cavity surface must belong to its enclosing metal body."
                )

        self.component_ids = np.empty(len(self.mesh.faces), dtype=int)
        for i, faces in enumerate(components):
            self.component_ids[faces] = i
