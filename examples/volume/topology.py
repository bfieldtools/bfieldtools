"""
Tree–cotree gauge and volume edge loops
=====================================

Load the bundled tetrahedral sphere and inspect independent current loops. Run
``uv run --locked --extra meshing python examples/volume/topology.py`` from the
repository root. Use the slider to select a retained edge. For an offscreen
preview, append ``--screenshot /tmp/volume-topology.png``.

The left panel distinguishes boundary edges (gray), interior tree edges
(orange, coefficients fixed to zero), and retained cotree edges (blue).
The boundary is contracted to a root when selecting the interior tree; the
orange edges need not form a connected tree before that contraction.

The right panel shows one retained edge (red), its incident tetrahedra, and
its constant cell currents for a coefficient of 1 A. Green segments form a
closed streamline through the cell centroids. Arrows show current direction
and relative magnitude, with a common display scale for the selected loop.
"""

import argparse

import numpy as np
import pyvista as pv

from bfieldtools.volume import TetraTopology, TreeCotreeGauge, load_example_tetra_mesh


def segments_polydata(segments):
    """VTK line segments from an array with shape (nsegments, 2, 3)."""
    segments = np.asarray(segments)
    points = segments.reshape(-1, 3)
    lines = np.column_stack(
        (np.full(len(segments), 2), np.arange(len(points)).reshape(-1, 2))
    )
    return pv.PolyData(points, lines=lines.ravel())


def make_plotter(mesh, topology, gauge, *, off_screen=False):
    """Create the two views; keep plotting separate from topology construction."""
    if not gauge.n_dofs:
        raise ValueError("the example needs at least one retained current loop")
    grid = mesh.to_pyvista()
    boundary = pv.PolyData(
        mesh.vertices,
        faces=np.column_stack(
            (np.full(len(mesh.boundary_faces), 3), mesh.boundary_faces)
        ).ravel(),
    )
    plotter = pv.Plotter(shape=(1, 2), window_size=(1400, 700), off_screen=off_screen)
    plotter.set_background("white")
    plotter.subplot(0, 0)
    plotter.add_text("Boundary, interior tree, and cotree", font_size=12)
    plotter.add_mesh(boundary, color="lightgray", opacity=0.06)
    for edge_ids, color, width, label in (
        (np.flatnonzero(topology.boundary_edge_mask), "gray", 1, "Boundary: zero"),
        (gauge.tree_edges, "darkorange", 5, "Interior tree: zero"),
        (gauge.retained_edges, "steelblue", 2, "Cotree: independent"),
    ):
        if len(edge_ids):
            plotter.add_mesh(
                segments_polydata(mesh.vertices[topology.edges[edge_ids]]),
                color=color,
                line_width=width,
                label=label,
            )
    plotter.add_legend(size=(0.48, 0.16), loc="lower left", face="line")
    plotter.view_isometric()
    plotter.add_axes()

    plotter.subplot(0, 1)
    plotter.add_mesh(boundary, color="lightgray", opacity=0.06)
    plotter.add_axes()

    def show_loop(value):
        selected = int(np.clip(round(value), 0, gauge.n_dofs - 1))
        edge_id = gauge.retained_edges[selected]
        edge = topology.edges[edge_id]
        cells, local_edges = np.nonzero(topology.cell_edges == edge_id)
        opposite = (local_edges + 3) % 6
        directions = mesh.cell_edge_vectors[cells, opposite]
        currents = (
            topology.cell_edge_signs[cells, local_edges, None]
            * directions
            / (3 * mesh.volumes[cells, None])
        )
        # lambda_i = lambda_j = 1/4 on this streamline. The other two
        # barycentric coordinates vary between (0, 1/2) and (1/2, 0).
        other_nodes = mesh.tetrahedra[cells]
        others = other_nodes[(other_nodes != edge[0]) & (other_nodes != edge[1])]
        others = others.reshape(-1, 2)
        loop = 0.5 * mesh.vertices[edge].mean(axis=0) + 0.5 * mesh.vertices[others]
        selected_line = segments_polydata(mesh.vertices[edge][None])

        plotter.subplot(0, 0)
        plotter.add_mesh(selected_line, color="crimson", line_width=7, name="selected")
        plotter.subplot(0, 1)
        plotter.add_text(
            f"Cotree loop {selected}: edge {edge[0]} → {edge[1]}\n"
            f"{len(cells)} incident cells; coefficient = 1 A",
            name="loop_title",
            font_size=12,
        )
        plotter.add_mesh(
            grid.extract_cells(cells),
            color="lightsteelblue",
            opacity=0.18,
            show_edges=True,
            name="support",
        )
        plotter.add_mesh(selected_line, color="crimson", line_width=7, name="edge")
        plotter.add_mesh(
            segments_polydata(loop),
            color="forestgreen",
            line_width=5,
            name="streamline",
        )
        plotter.add_arrows(
            mesh.cell_centers[cells],
            currents,
            mag=0.25 / np.linalg.norm(currents, axis=1).max(),
            color="black",
            name="currents",
        )

    midpoints = mesh.vertices[topology.edges[gauge.retained_edges]].mean(axis=1)
    selected = np.argmin(
        np.linalg.norm(midpoints - mesh.cell_centers.mean(axis=0), axis=1)
    )
    show_loop(selected)
    plotter.view_isometric()
    plotter.link_views()
    if gauge.n_dofs > 1:
        plotter.add_slider_widget(
            show_loop,
            rng=(0, gauge.n_dofs - 1),
            value=int(selected),
            title="Retained edge index",
            fmt="%.0f",
            pointa=(0.1, 0.08),
            pointb=(0.9, 0.08),
            interaction_event="end",
        )
    return plotter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screenshot", help="render off screen to this image path")
    args = parser.parse_args()
    mesh = load_example_tetra_mesh()
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    print(f"Cells: {len(mesh.tetrahedra)}; edges: {len(topology.edges)}")
    print(f"Boundary edges: {np.count_nonzero(topology.boundary_edge_mask)}")
    print(
        f"Interior tree edges: {len(gauge.tree_edges)}; independent loops: {gauge.n_dofs}"
    )
    for name, residual in (
        ("C G", topology.curl @ topology.gradient),
        ("D C", topology.divergence @ topology.curl),
        ("boundary flux", gauge.curl[topology.boundary_face_mask]),
    ):
        residual.eliminate_zeros()
        print(f"{name}: {residual.nnz} nonzero residuals")
        if residual.nnz:
            raise RuntimeError(f"failed topology identity: {name}")
    plotter = make_plotter(mesh, topology, gauge, off_screen=bool(args.screenshot))
    plotter.show(screenshot=args.screenshot)


if __name__ == "__main__":
    main()
