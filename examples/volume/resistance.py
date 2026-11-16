"""
Resistance spectrum and volume-current eigenfunctions
====================================================

Run ``uv run --locked --extra meshing python examples/volume/resistance.py``
from the repository root. The figure shows the spectrum of R and three
eigenfunctions on a tetrahedral sphere. Use ``--modes 0 5 20`` to choose modes
(zero-based indices). Append ``--screenshot /tmp/volume-resistance.png`` for an
offscreen image. ``--resistivity`` sets a uniform scalar value in ohm metres;
the default is 1, and the sphere radius is 1 metre.

These are ordinary coefficient-space eigenpairs R v = lambda v, with Euclidean
normalization ||v||_2 = 1 A. They depend on the selected tree gauge. Their
eigenvalues measure Joule loss per squared coefficient norm; they are not decay
rates. Physical decay modes will require the inductance metric R v = lambda M v.

Colors show |J| on a central slice. Black arrows show the constant cell
current directions in the front half of the sphere, with fixed display lengths.
The three panels share a color scale. The translucent
outer surface gives context. Current tangential components may jump between
cells, while normal fluxes match across faces.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from scipy.linalg import eigh
import trimesh

from bfieldtools.volume import (
    TreeCotreeGauge,
    build_edge_current_basis,
    TetraTopology,
    reconstruct_currents,
    resistance_matrix,
    tetrahedralize,
)


def make_plotter(basis, eigenvalues, eigenvectors, modes, *, off_screen=False):
    """Show a spectrum chart and selected current eigenfunctions in one window."""
    mesh = basis.mesh
    fields = reconstruct_currents(basis, eigenvectors[:, modes])
    magnitudes = np.linalg.norm(fields, axis=1)
    max_current = magnitudes.max()
    colors = ("tab:blue", "tab:orange", "tab:green")
    figure, axes = plt.subplots(figsize=(6, 5))
    axes.semilogy(np.arange(len(eigenvalues)), eigenvalues, color="0.35", linewidth=1)
    for mode, color in zip(modes, colors):
        axes.scatter(
            mode, eigenvalues[mode], color=color, label=f"Mode {mode}", zorder=3
        )
    axes.set(
        xlabel="Mode index (ascending)",
        ylabel="Eigenvalue (Ω)",
        title="Resistance: R v = λ v, ‖v‖₂ = 1 A",
    )
    axes.grid(True, which="both", alpha=0.2)
    axes.legend()
    figure.text(
        0.5,
        0.01,
        "Spectrum depends on tree gauge; decay modes require M.",
        ha="center",
        fontsize=9,
    )
    figure.subplots_adjust(left=0.14, bottom=0.18, right=0.97, top=0.91)

    plotter = pv.Plotter(shape=(2, 2), window_size=(1300, 1000), off_screen=off_screen)
    plotter.set_background("white")
    plotter.subplot(0, 0)
    plotter.add_chart(pv.ChartMPL(figure))
    # The VTK chart owns the figure; suppress a second GUI window.
    plt.close(figure)
    boundary = pv.PolyData(
        mesh.vertices,
        faces=np.column_stack(
            (np.full(len(mesh.boundary_faces), 3), mesh.boundary_faces)
        ).ravel(),
    )
    grid = mesh.to_pyvista()
    center = mesh.cell_centers.mean(axis=0)
    front = mesh.cell_centers[:, 0] >= center[0]
    radius = np.ptp(mesh.vertices, axis=0).max() / 2
    for panel, mode in enumerate(modes):
        subplot = panel + 1
        plotter.subplot(subplot // 2, subplot % 2)
        plotter.add_text(
            f"Mode {mode}: eigenvalue = {eigenvalues[mode]:.4g} ohm",
            font_size=12,
            color=colors[panel],
        )
        plotter.add_mesh(boundary, color="lightgray", opacity=0.04)
        grid.cell_data["current_magnitude"] = magnitudes[:, panel]
        plotter.add_mesh(
            grid.slice(normal="x", origin=center),
            scalars="current_magnitude",
            clim=(0, max_current),
            cmap="YlOrRd",
            show_edges=True,
            scalar_bar_args={"title": "|J| (A/m^2)", "position_x": 0.12, "width": 0.76},
        )
        vectors = fields[front, :, panel]
        lengths = np.linalg.norm(vectors, axis=1)
        directions = np.divide(
            vectors,
            lengths[:, None],
            out=np.zeros_like(vectors),
            where=lengths[:, None] > 0,
        )
        plotter.add_arrows(
            mesh.cell_centers[front], directions, mag=0.15 * radius, color="black"
        )
        plotter.add_axes()
        plotter.view_vector((3, 1, 1))
    plotter.link_views(views=(1, 2, 3))
    return plotter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resistivity", type=float, default=1.0)
    parser.add_argument("--modes", nargs=3, type=int, default=[0, 1, 2])
    parser.add_argument("--screenshot", help="render off screen to this image path")
    args = parser.parse_args()
    mesh = tetrahedralize(trimesh.creation.icosphere(subdivisions=1), max_volume=0.05)
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    basis = build_edge_current_basis(mesh, topology, gauge)
    matrix = resistance_matrix(basis, args.resistivity)
    if min(args.modes) < 0 or max(args.modes) >= basis.n_dofs:
        parser.error(f"mode indices must be between 0 and {basis.n_dofs - 1}")
    eigenvalues, eigenvectors = eigh(matrix.toarray())
    if eigenvalues[0] <= 0:
        raise RuntimeError("resistance is not positive definite")
    residual = np.linalg.norm(
        matrix @ eigenvectors - eigenvectors * eigenvalues
    ) / np.linalg.norm(matrix.toarray())
    currents = reconstruct_currents(basis, eigenvectors[:, args.modes])
    losses = args.resistivity * np.sum(
        mesh.volumes[:, None] * np.sum(currents**2, axis=1), axis=0
    )
    np.testing.assert_allclose(losses, eigenvalues[args.modes], rtol=1e-10)
    print(f"Cells: {len(mesh.tetrahedra)}; independent coefficients: {basis.n_dofs}")
    print(f"Resistance: {matrix.shape}, {matrix.nnz} nonzeros")
    print(f"Eigenvalues: {eigenvalues[0]:.6g} … {eigenvalues[-1]:.6g} Ω")
    print(f"Relative eigenpair residual: {residual:.3g}")
    print(f"Selected-mode Joule losses (W): {losses}")
    plotter = make_plotter(
        basis, eigenvalues, eigenvectors, args.modes, off_screen=bool(args.screenshot)
    )
    plotter.show(screenshot=args.screenshot)


if __name__ == "__main__":
    main()
