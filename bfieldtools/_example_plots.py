"""PyVista dataset constructors shared by the example scripts.

These helpers only prepare geometry; examples own their Plotter and actors.
"""

import numpy as np
import pyvista as pv


def point_cloud(x, y, z):
    """Create points from coordinate arrays."""
    return pv.PolyData(
        np.column_stack((np.ravel(x), np.ravel(y), np.ravel(z))).astype(float)
    )


def triangle_mesh(x, y, z, triangles):
    """Create a triangular surface from coordinates and face indices."""
    faces = np.column_stack((np.full(len(triangles), 3), triangles))
    return pv.PolyData(point_cloud(x, y, z).points, faces)


def polyline(x, y, z, scalars=None):
    """Create one connected line, optionally coloured by a point scalar."""
    points = point_cloud(x, y, z)
    points.lines = np.r_[points.n_points, np.arange(points.n_points)]
    if scalars is not None:
        points["scalars"] = np.ravel(scalars)
    return points


def vector_glyphs(
    x, y, z, u, v, w, scale_factor=None, mode="arrow", scalars=None, scale_mode="vector"
):
    """Build oriented glyphs, preserving vector magnitude and optional colours."""
    points = point_cloud(x, y, z)
    points["vectors"] = np.column_stack((np.ravel(u), np.ravel(v), np.ravel(w)))
    points["magnitude"] = np.linalg.norm(points["vectors"], axis=1)
    if scalars is not None:
        points["scalars"] = np.ravel(scalars)
    if scale_factor is None:
        maximum = points["magnitude"].max(initial=0)
        scale_factor = (points.length or 1) * 0.05 / maximum if maximum else 1
    geom = pv.Arrow()
    if mode == "cylinder":
        geom = pv.Cylinder(direction=(1, 0, 0), radius=0.1, height=1)
    elif mode in ("2ddash", "2darrow"):
        geom = pv.Line((-0.5, 0, 0), (0.5, 0, 0))
    scale = False if scale_mode == "none" else "magnitude"
    if scale_mode == "scalar" and scalars is not None:
        scale = "scalars"
    glyphs = points.glyph(orient="vectors", scale=scale, factor=scale_factor, geom=geom)
    glyphs.set_active_scalars("scalars" if scalars is not None else "magnitude")
    return glyphs


def structured_field(x, y, z, u, v, w, scalars=None):
    """Build a structured grid using VTK's Fortran point ordering."""
    grid = pv.StructuredGrid(x, y, z)
    grid["vectors"] = np.column_stack(
        [np.asarray(a).ravel(order="F") for a in (u, v, w)]
    )
    grid["magnitude"] = np.linalg.norm(grid["vectors"], axis=1)
    if scalars is not None:
        grid["scalars"] = np.asarray(scalars).ravel(order="F")
    grid.set_active_vectors("vectors")
    grid.set_active_scalars("magnitude")
    return grid


def set_view(plotter, azimuth=0, elevation=0, distance=None, roll=None):
    """Position the camera using spherical angles in degrees."""
    center = np.asarray(plotter.center)
    if distance is None:
        distance = max(plotter.length * 2, 1)
    azimuth, elevation = np.deg2rad([azimuth, elevation])
    offset = distance * np.array(
        [
            np.sin(elevation) * np.cos(azimuth),
            np.sin(elevation) * np.sin(azimuth),
            np.cos(elevation),
        ]
    )
    plotter.camera_position = (
        center + offset,
        center,
        (0, 1, 0) if abs(np.sin(elevation)) < 1e-8 else (0, 0, 1),
    )
    if roll is not None:
        plotter.camera.roll = roll


def mesh_actor(figure):
    """Return the last mesh actor from a Plotter, or an already supplied Actor."""
    if isinstance(figure, pv.Actor):
        return figure
    return next(
        actor
        for actor in reversed(list(figure.actors.values()))
        if isinstance(actor, pv.Actor)
    )
