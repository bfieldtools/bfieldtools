"""Backend-independent linear tetrahedral meshes and optional TetGen meshing.

Coordinates retain the input length unit; volumes use its cube. No geometry
repair, rescaling, or automatic changes to surface topology are performed.
"""

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray
import trimesh

if TYPE_CHECKING:
    import pyvista as pv


# Positive tetrahedron: each edge curl points along edge (n + 3) % 6.
TETRA_EDGES = ((0, 1), (0, 2), (0, 3), (2, 3), (3, 1), (1, 2))
TETRA_FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))


def _indices(
    values: ArrayLike, width: int, nvertices: int, name: str
) -> NDArray[np.int64]:
    values = np.asarray(values)
    if values.ndim != 2 or values.shape[1] != width or not len(values):
        raise ValueError(f"{name} must have nonempty shape (n, {width})")
    if not np.issubdtype(values.dtype, np.integer):
        raise ValueError(f"{name} must contain integer vertex indices")
    if values.min() < 0 or values.max() >= nvertices:
        raise ValueError(f"{name} contains an out-of-range vertex index")
    values = np.array(values, dtype=np.int64, copy=True)
    if np.any(np.diff(np.sort(values, axis=1), axis=1) == 0):
        raise ValueError(f"{name} contains repeated vertices")
    if len(np.unique(np.sort(values, axis=1), axis=0)) != len(values):
        raise ValueError(f"{name} contains duplicate elements")
    return values


def _vertices(values: ArrayLike) -> NDArray[np.float64]:
    values = np.array(values, dtype=float, copy=True)
    if values.ndim != 2 or values.shape[1] != 3 or not len(values):
        raise ValueError("vertices must have nonempty shape (n, 3)")
    if not np.isfinite(values).all():
        raise ValueError("vertices must be finite")
    return values


class TetraMesh:
    """A validated, positively oriented mesh of four-node tetrahedra.

    Parameters
    ----------
    vertices : array-like, shape (nvertices, 3)
        Vertex coordinates.
    tetrahedra : integer array-like, shape (ncells, 4)
        Zero-based vertex indices. Negative orientations are corrected on a
        copy. Degenerate cells, duplicate cells, and nonconforming shared-face
        orientations are rejected.

    Notes
    -----
    Arrays are copied and made read-only. Validation checks local geometry and
    face incidence, not arbitrary geometric overlaps or hanging nodes. Supplied
    meshes must already be conforming. No meshing backend is needed to construct
    this container.
    """

    vertices: NDArray[np.float64]
    tetrahedra: NDArray[np.int64]
    volumes: NDArray[np.float64]
    boundary_faces: NDArray[np.int64]

    def __init__(self, vertices: ArrayLike, tetrahedra: ArrayLike) -> None:
        self.vertices = _vertices(vertices)
        self.tetrahedra = _indices(tetrahedra, 4, len(self.vertices), "tetrahedra")
        points = self.vertices[self.tetrahedra]
        edges = points[:, 1:] - points[:, :1]
        determinants = np.linalg.det(edges)
        scale = np.max(np.linalg.norm(edges, axis=2), axis=1)
        if np.any(np.abs(determinants) <= 64 * np.finfo(float).eps * scale**3):
            raise ValueError("tetrahedra contains degenerate or near-degenerate cells")
        negative = determinants < 0
        self.tetrahedra[negative] = self.tetrahedra[negative][:, [1, 0, 2, 3]]
        self.volumes = np.abs(determinants) / 6

        # For positive determinants these local face orderings point outward.
        faces = self.tetrahedra[:, TETRA_FACES]
        faces = faces.reshape(-1, 3)
        _, first, inverse, counts = np.unique(
            np.sort(faces, axis=1),
            axis=0,
            return_index=True,
            return_inverse=True,
            return_counts=True,
        )
        if np.any(counts > 2):
            raise ValueError("nonmanifold mesh: a face belongs to more than two cells")
        # Shared faces must have opposite orientation (permutation parity).
        parity = np.prod(np.sign(faces[:, [1, 2, 2]] - faces[:, [0, 0, 1]]), axis=1)
        balance = np.bincount(inverse, weights=parity)
        if np.any(balance[counts == 2] != 0):
            raise ValueError("shared faces must have opposite outward orientations")
        self.boundary_faces = faces[first[counts == 1]]
        for array in (
            self.vertices,
            self.tetrahedra,
            self.volumes,
            self.boundary_faces,
        ):
            array.flags.writeable = False

    @cached_property
    def cell_centers(self) -> NDArray[np.float64]:
        """Cell centroids, shape (ncells, 3)."""
        centers = self.vertices[self.tetrahedra].mean(axis=1)
        centers.flags.writeable = False
        return centers

    @cached_property
    def cell_edge_vectors(self) -> NDArray[np.float64]:
        """Directed edge vectors, shape (ncells, 6, 3), in ``TETRA_EDGES`` order.

        The constant local edge curls are the opposing vectors divided by
        ``3 * volumes``. Local edge n has opposing edge ``(n + 3) % 6``.
        """
        edges = np.asarray(TETRA_EDGES)
        vectors = (
            self.vertices[self.tetrahedra[:, edges[:, 1]]]
            - self.vertices[self.tetrahedra[:, edges[:, 0]]]
        )
        vectors.flags.writeable = False
        return vectors

    @cached_property
    def cell_quality(self) -> NDArray[np.float64]:
        """Mean-ratio quality: 1 for a regular tetrahedron, tending to 0 for slivers."""
        points = self.vertices[self.tetrahedra]
        edge_sum = sum(
            np.sum((points[:, i] - points[:, j]) ** 2, axis=1)
            for i in range(4)
            for j in range(i)
        )
        quality = 12 * (3 * self.volumes) ** (2 / 3) / edge_sum
        quality.flags.writeable = False
        return quality

    def to_pyvista(self) -> pv.UnstructuredGrid:
        """Return a PyVista UnstructuredGrid with volume and quality cell data."""
        import pyvista as pv

        cells = np.column_stack((np.full(len(self.tetrahedra), 4), self.tetrahedra))
        grid = pv.UnstructuredGrid(
            cells.ravel(),
            np.full(len(cells), pv.CellType.TETRA, dtype=np.uint8),
            self.vertices.copy(),
        )
        grid.cell_data["volume"] = self.volumes.copy()
        grid.cell_data["quality"] = self.cell_quality
        return grid


def tetrahedralize(
    surface: trimesh.Trimesh,
    *,
    max_volume: float | None = None,
    min_ratio: float = 2.0,
    min_dihedral: float = 0.0,
) -> TetraMesh:
    """Fill one closed, consistently oriented triangular surface using TetGen.

    Parameters
    ----------
    surface : trimesh.Trimesh
        One outward-oriented, watertight genus-zero boundary. STL input should
        be loaded with trimesh's vertex merging enabled (the default). Cavities,
        multiple solids, material interfaces, and boundary labels are not yet
        supported. Self-intersections must be excluded by the caller; the checks
        here establish topology and local validity, not geometric intersection.
    max_volume : float, optional
        Requested maximum tetrahedron volume, in cubed input length units.
        Enables TetGen's fixed-volume refinement. Smaller values produce finer
        meshes, but TetGen can return cells exceeding this target; inspect
        ``mesh.volumes`` when a strict bound is required.
    min_ratio : float, default 2.0
        Maximum circumradius / shortest-edge ratio. Must exceed 1; smaller values
        request better quality and can greatly increase meshing cost.
    min_dihedral : float, default 0
        Requested minimum dihedral angle in degrees, in [0, 60). Aggressive
        quality constraints can prevent the mesher from terminating.

    Returns
    -------
    TetraMesh
        First-order cells with positive orientation and outward boundary faces.

    Notes
    -----
    Install the optional backend with ``pip install 'bfieldtools[meshing]'``.
    TetGen may insert boundary and interior vertices, preserving the piecewise
    planar boundary rather than projecting onto an ideal curved surface.
    """
    if not isinstance(surface, trimesh.Trimesh):
        raise TypeError("surface must be a trimesh.Trimesh")
    vertices = _vertices(surface.vertices)
    faces = _indices(surface.faces, 3, len(vertices), "surface faces")
    triangles = vertices[faces]
    edges = triangles[:, 1:] - triangles[:, :1]
    areas2 = np.linalg.norm(np.cross(edges[:, 0], edges[:, 1]), axis=1)
    scale = np.max(np.linalg.norm(edges, axis=2), axis=1)
    if np.any(areas2 <= 64 * np.finfo(float).eps * scale**2):
        raise ValueError("surface contains degenerate triangles")
    checked = trimesh.Trimesh(vertices, faces, process=False)
    if not checked.is_watertight:
        raise ValueError("surface must be watertight with two triangles per edge")
    if not checked.is_winding_consistent:
        raise ValueError("surface must have consistent triangle winding")
    if len(checked.split(only_watertight=False)) != 1 or checked.euler_number != 2:
        raise ValueError("only a single genus-zero boundary is currently supported")
    if checked.volume <= 0:
        raise ValueError("surface must face outward; invert its winding explicitly")
    if max_volume is not None and (not np.isfinite(max_volume) or max_volume <= 0):
        raise ValueError("max_volume must be finite and positive")
    if not np.isfinite(min_ratio) or min_ratio <= 1:
        raise ValueError("min_ratio must be finite and greater than 1")
    if not np.isfinite(min_dihedral) or not 0 <= min_dihedral < 60:
        raise ValueError("min_dihedral must be finite and in [0, 60)")
    try:
        import tetgen
    except ImportError as exc:
        raise ImportError(
            "Tetrahedralization requires TetGen: pip install 'bfieldtools[meshing]'"
        ) from exc

    backend = tetgen.TetGen(vertices, faces.astype(np.int32))
    result = backend.tetrahedralize(
        plc=True,
        order=1,
        quality=True,
        minratio=min_ratio,
        mindihedral=min_dihedral,
        fixedvolume=max_volume is not None,
        maxvolume=-1.0 if max_volume is None else max_volume,
    )
    mesh = TetraMesh(result[0], result[1])
    if not np.isclose(mesh.volumes.sum(), checked.volume, rtol=1e-7, atol=0):
        raise ValueError("tetrahedral mesh volume does not match the input surface")
    return mesh
