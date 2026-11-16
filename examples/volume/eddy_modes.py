"""
Inductance and freely decaying volume-current modes
=================================================

Run ``uv run --locked --extra meshing python examples/volume/eddy_modes.py``
from the repository root. The figure shows decay rates and the first three
current eigenfunctions of ``R v = lambda M v`` on the bundled 808-cell tetrahedral sphere.
``--modes 0 3 8`` chooses three zero-based indices. Smaller rates correspond to
longer decay times ``tau = 1/lambda``.

The default sphere radius is 0.1 m, with uniform resistivity 1e-8 ohm metres.
Both can be changed with ``--radius`` and ``--resistivity``. Geometry is loaded
from ``example_meshes/unit_sphere_tetra.npz``; TetGen is not needed at runtime.
``--quadrature-rtol`` controls successive outer quadrature refinements, not
spatial discretization error. ``--check-convergence`` repeats assembly at a
tighter integration tolerance and reports changes in the selected decay rates.

Eigenvectors are M-orthonormal. For display, each current magnitude is divided
by its own peak. Colors show this relative magnitude on a central slice, and
fixed-length black arrows show current directions in the front half of the
sphere. Rotate the three linked views to inspect their spatial patterns.
Modes in a nearly degenerate group can mix; their span is more meaningful than
an individual orientation. The mesh remains a demonstration rather than a spatial convergence study.

Append ``--screenshot /tmp/volume-eddy-modes.png`` for an offscreen image.
"""

import argparse
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from scipy.linalg import eigh

from bfieldtools.volume import (
    TetraTopology,
    TreeCotreeGauge,
    build_edge_current_basis,
    inductance_matrix,
    reconstruct_currents,
    resistance_matrix,
    TetraMesh,
    load_example_tetra_mesh,
)


def make_plotter(basis, rates, vectors, modes, *, off_screen=False):
    """Create the spectrum and linked views of three current eigenfunctions."""
    mesh = basis.mesh
    fields = reconstruct_currents(basis, vectors[:, modes])
    magnitudes = np.linalg.norm(fields, axis=1)
    relative = magnitudes / magnitudes.max(axis=0)
    colors = ("tab:blue", "tab:orange", "tab:green")
    figure, axes = plt.subplots(figsize=(6, 5))
    axes.semilogy(np.arange(len(rates)), rates, color="0.35", linewidth=1)
    for mode, color in zip(modes, colors):
        axes.scatter(mode, rates[mode], color=color, label=f"Mode {mode}", zorder=3)
    axes.set(
        xlabel="Mode index (ascending rate)",
        ylabel="Decay rate λ (s⁻¹)",
        title="Volume eddy currents: R v = λ M v",
    )
    axes.grid(True, which="both", alpha=0.2)
    axes.legend()
    figure.text(
        0.5,
        0.015,
        "Lowest rates give the longest-lived currents.",
        ha="center",
        fontsize=9,
    )
    figure.subplots_adjust(left=0.15, bottom=0.18, right=0.97, top=0.91)

    plotter = pv.Plotter(shape=(2, 2), window_size=(1300, 1000), off_screen=off_screen)
    plotter.set_background("white")
    plotter.subplot(0, 0)
    plotter.add_chart(pv.ChartMPL(figure))
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
            f"Mode {mode}: rate = {rates[mode]:.4g} /s\n"
            f"Decay time = {1e3/rates[mode]:.4g} ms",
            font_size=12,
            color=colors[panel],
        )
        plotter.add_mesh(boundary, color="lightgray", opacity=0.04)
        grid.cell_data["relative_magnitude"] = relative[:, panel]
        plotter.add_mesh(
            grid.slice(normal="x", origin=center),
            scalars="relative_magnitude",
            clim=(0, 1),
            cmap="YlOrRd",
            show_edges=True,
            lighting=False,
            scalar_bar_args={
                "title": "|J| / mode peak",
                "position_x": 0.14,
                "width": 0.72,
            },
        )
        current = fields[front, :, panel]
        lengths = np.linalg.norm(current, axis=1)
        directions = np.divide(
            current,
            lengths[:, None],
            out=np.zeros_like(current),
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
    parser.add_argument("--radius", type=float, default=0.1)
    parser.add_argument("--resistivity", type=float, default=1e-8)
    parser.add_argument("--quadrature-rtol", type=float, default=1e-4)
    parser.add_argument("--max-order", type=int, default=12)
    parser.add_argument("--check-convergence", action="store_true")
    parser.add_argument("--modes", type=int, nargs=3, default=[0, 1, 2])
    parser.add_argument("--screenshot", help="render off screen to this image path")
    args = parser.parse_args()
    if not np.isfinite(args.radius) or args.radius <= 0:
        parser.error("radius must be positive and finite")
    unit = load_example_tetra_mesh()
    mesh = TetraMesh(unit.vertices * args.radius, unit.tetrahedra)
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    basis = build_edge_current_basis(mesh, topology, gauge)
    if min(args.modes) < 0 or max(args.modes) >= basis.n_dofs:
        parser.error(f"mode indices must be between 0 and {basis.n_dofs-1}")
    resistance = resistance_matrix(basis, args.resistivity).toarray()
    print(
        f"Assembling inductance for {len(mesh.tetrahedra)} cells; this may take a few minutes.",
        flush=True,
    )
    start = perf_counter()
    magnetic = inductance_matrix(
        basis, rtol=args.quadrature_rtol, max_order=args.max_order
    )
    print(f"Cells: {len(mesh.tetrahedra)}; independent currents: {basis.n_dofs}")
    print(f"Inductance assembly: {perf_counter()-start:.2f} s")
    rates, vectors = eigh(resistance, magnetic)
    if rates[0] <= 0:
        raise RuntimeError("decay rates must be positive")
    residual = resistance @ vectors - (magnetic @ vectors) * rates
    relative_residual = np.linalg.norm(residual) / np.linalg.norm(resistance @ vectors)
    orthogonality = np.linalg.norm(vectors.T @ magnetic @ vectors - np.eye(len(rates)))
    print(f"Relative eigenpair residual: {relative_residual:.3g}")
    print(f"M-orthogonality error: {orthogonality:.3g}")
    if relative_residual > 1e-10 or orthogonality > 1e-10:
        raise RuntimeError("generalized eigenpair verification failed")
    for mode in args.modes:
        print(f"Mode {mode}: rate={rates[mode]:.6g} /s, tau={1e3/rates[mode]:.6g} ms")
    if args.check_convergence:
        refined = inductance_matrix(
            basis, rtol=args.quadrature_rtol / 4, max_order=args.max_order
        )
        refined_rates = eigh(
            resistance, refined, eigvals_only=True, subset_by_index=(0, max(args.modes))
        )
        shifts = (
            np.abs(refined_rates[args.modes] - rates[args.modes])
            / refined_rates[args.modes]
        )
        print(f"Selected-rate changes after tightening quadrature: {shifts}")
    plotter = make_plotter(
        basis, rates, vectors, args.modes, off_screen=bool(args.screenshot)
    )
    plotter.show(screenshot=args.screenshot)


if __name__ == "__main__":
    main()
