"""
Tetrahedral mesh of a solid sphere
================================

Fill the bundled unit_sphere.stl boundary with linear tetrahedra and inspect
their interiors. Install with ``pip install 'bfieldtools[meshing]'``, or run
``uv run --extra meshing python examples/volume/sphere_meshing.py`` from the
repository root. STL loading merges coincident vertices; no holes are filled
and no boundary geometry is repaired by the meshing utility.

The mesh fills the faceted STL, so its volume differs slightly from 4*pi/3.
The cutaway removes whole cells to expose tetrahedral faces. The second panel
shows a planar cross section colored by mean-ratio cell quality (1 is ideal).
"""

import numpy as np
import pyvista as pv

from bfieldtools.utils import load_example_mesh
from bfieldtools.volume import tetrahedralize


surface = load_example_mesh("unit_sphere")
mesh = tetrahedralize(surface, max_volume=0.01)
grid = mesh.to_pyvista()

print(f"Vertices: {len(mesh.vertices)}, tetrahedra: {len(mesh.tetrahedra)}")
print(f"Boundary triangles: {len(mesh.boundary_faces)}")
print(f"Mesh volume: {mesh.volumes.sum():.6f}; STL volume: {surface.volume:.6f}")
print(f"Ideal sphere volume: {4 * np.pi / 3:.6f}")
print(f"Cell volume range: {mesh.volumes.min():.3g} – {mesh.volumes.max():.3g}")
print(f"Minimum mean-ratio quality: {mesh.cell_quality.min():.3f}")

# Select complete cells, exposing the triangular faces inside the solid.
center = surface.bounds.mean(axis=0)
cutaway = grid.extract_cells(mesh.cell_centers[:, 0] < center[0])
boundary = pv.PolyData(
    surface.vertices, np.column_stack((np.full(len(surface.faces), 3), surface.faces))
)

plotter = pv.Plotter(shape=(1, 2), window_size=(1200, 600))
plotter.subplot(0, 0)
plotter.add_text("Solid sphere: tetrahedral cutaway", font_size=12)
plotter.add_mesh(boundary, color="lightgray", opacity=0.12)
plotter.add_mesh(cutaway, scalars="quality", clim=(0, 1), show_edges=True)
plotter.view_isometric()
plotter.add_axes()

plotter.subplot(0, 1)
plotter.add_text("Central slice: cell quality", font_size=12)
plotter.add_mesh(
    grid.slice(normal="x", origin=center),
    scalars="quality",
    clim=(0, 1),
    show_edges=True,
)
plotter.view_yz()
plotter.add_axes()
plotter.show()
