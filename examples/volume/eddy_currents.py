"""Eddy currents in a solid sphere driven by a nearby prescribed current loop.

Run from the repository root::

    uv run --locked --extra meshing python examples/volume/eddy_currents.py

The bundled 808-cell unit sphere is scaled to radius 0.1 m. Resistivity is
1e-8 ohm metres, and a 1000 A smooth pulse flows in a 0.075 m radius polygonal
loop at z=0.14 m. All quantities use SI units. The conductor starts unexcited;
after the 80 ms pulse, its currents decay freely. All generalized modes are
used, so there is no modal truncation. The loop current is prescribed.

The figure shows the source waveform, decay times, modal amplitudes, and three
current snapshots with a common physical magnitude scale. Black arrows show
directions at fixed display length. ``--check-convergence`` tightens inductance
quadrature and refines coupling integration and loop segmentation. These are
integration checks; spatial convergence still requires a mesh-refinement study.
``--screenshot /tmp/volume-eddy-currents.png`` renders off screen.
"""

import argparse
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from scipy.integrate import solve_ivp
from scipy.linalg import cho_factor, cho_solve, eigh

from bfieldtools.volume import (
    TetraMesh,
    TetraTopology,
    TreeCotreeGauge,
    build_edge_current_basis,
    inductance_matrix,
    load_example_tetra_mesh,
    loop_coupling,
    modal_response,
    reconstruct_currents,
    resistance_matrix,
)
from bfieldtools.volume.mesh import TETRA_FACES


def source_loop(segments=96):
    angle = np.linspace(0, 2 * np.pi, segments + 1)
    vertices = np.column_stack(
        (0.075 * np.cos(angle), 0.075 * np.sin(angle), np.full_like(angle, 0.14))
    )
    vertices[-1] = vertices[0]
    return vertices


PULSE_DURATION = 0.08
PEAK_CURRENT = 1000.0


def waveform(t):
    t = np.asarray(t)
    return np.where(
        (t > 0) & (t < PULSE_DURATION),
        PEAK_CURRENT * np.sin(np.pi * t / PULSE_DURATION) ** 2,
        0.0,
    )


def current_derivative(t):
    if 0 < t < PULSE_DURATION:
        return (
            PEAK_CURRENT
            * np.pi
            / PULSE_DURATION
            * np.sin(2 * np.pi * t / PULSE_DURATION)
        )
    return 0.0


def make_plotter(
    basis, loop, rates, times, amplitudes, fields, snapshots, *, off_screen=False
):
    """Plot the drive, modal dynamics, and three physical current snapshots."""
    plotter = pv.Plotter(shape=(2, 3), window_size=(1600, 1000), off_screen=off_screen)
    plotter.set_background("white")
    selected = np.argsort(np.max(np.abs(amplitudes), axis=1))[-4:][::-1]
    for panel in range(3):
        figure, axes = plt.subplots(figsize=(5, 4))
        if panel == 0:
            axes.plot(1e3 * times, waveform(times), color="tab:red")
            for t in snapshots:
                axes.axvline(1e3 * t, color="0.6", linestyle=":")
            axes.set(
                xlabel="Time (ms)",
                ylabel="Loop current (A)",
                title="Prescribed smooth current pulse",
            )
        elif panel == 1:
            axes.semilogy(np.arange(len(rates)), 1e3 / rates, color="tab:blue")
            axes.set(
                xlabel="Mode index",
                ylabel="Decay time (ms)",
                title="R v = λ M v; τ = 1/λ",
            )
        else:
            for mode in selected:
                axes.plot(1e3 * times, amplitudes[mode], label=f"Mode {mode}")
            axes.axvline(1e3 * PULSE_DURATION, color="0.6", linestyle=":")
            axes.set(
                xlabel="Time (ms)",
                ylabel="Modal amplitude (√J)",
                title="Four most strongly excited modes",
            )
            axes.legend()
        axes.grid(alpha=0.2)
        figure.tight_layout()
        plotter.subplot(0, panel)
        plotter.add_chart(pv.ChartMPL(figure))
        plt.close(figure)

    mesh = basis.mesh
    grid = mesh.to_pyvista()
    boundary = grid.extract_surface(algorithm="dataset_surface")
    magnitudes = np.linalg.norm(fields, axis=1)
    maximum = magnitudes.max()
    front = mesh.cell_centers[:, 0] >= 0
    wire = pv.lines_from_points(loop)
    for panel, t in enumerate(snapshots):
        plotter.subplot(1, panel)
        plotter.add_text(
            f"t = {1e3*t:.0f} ms; I = {float(waveform(t)):.0f} A", font_size=12
        )
        plotter.add_mesh(boundary, color="lightgray", opacity=0.04)
        plotter.add_mesh(wire, color="tab:red", line_width=5)
        grid.cell_data["current"] = magnitudes[:, panel]
        plotter.add_mesh(
            grid.slice(normal="x", origin=(0, 0, 0)),
            scalars="current",
            clim=(0, maximum),
            cmap="YlOrRd",
            lighting=False,
            show_edges=True,
            scalar_bar_args={"title": "|J| (A/m²)", "position_x": 0.15, "width": 0.7},
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
            mesh.cell_centers[front], directions, mag=0.012, color="black"
        )
        plotter.add_axes()
        plotter.view_vector((3, 1, 1))
    plotter.link_views(views=(3, 4, 5))
    return plotter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-convergence", action="store_true")
    parser.add_argument("--screenshot")
    args = parser.parse_args()
    unit = load_example_tetra_mesh()
    mesh = TetraMesh(0.1 * unit.vertices, unit.tetrahedra)
    topology = TetraTopology.from_mesh(mesh)
    gauge = TreeCotreeGauge.from_topology(topology)
    basis = build_edge_current_basis(mesh, topology, gauge)
    resistance = resistance_matrix(basis, 1e-8).toarray()
    print(
        f"Assembling inductance for {len(mesh.tetrahedra)} cells; this may take a few minutes.",
        flush=True,
    )
    start = perf_counter()
    magnetic = inductance_matrix(basis)
    rates, vectors = eigh(resistance, magnetic)
    print(
        f"Cells: {len(mesh.tetrahedra)}; currents: {basis.n_dofs}; assembly and modes: {perf_counter()-start:.1f} s",
        flush=True,
    )
    loop = source_loop()
    coupling = loop_coupling(basis, loop)
    refined_coupling = loop_coupling(basis, loop, quadrature_order=6)
    coupling_error = np.linalg.norm(
        vectors.T @ (refined_coupling - coupling)
    ) / np.linalg.norm(vectors.T @ refined_coupling)
    print(f"Coupling quadrature relative change: {coupling_error:.3g}")
    if coupling_error > 1e-4:
        raise RuntimeError("coupling quadrature needs further refinement")
    residual = np.linalg.norm(
        resistance @ vectors - (magnetic @ vectors) * rates
    ) / np.linalg.norm(resistance @ vectors)
    orthogonality = np.linalg.norm(vectors.T @ magnetic @ vectors - np.eye(len(rates)))
    if rates[0] <= 0 or residual > 1e-10 or orthogonality > 1e-9:
        raise RuntimeError("generalized modes failed verification")
    print(f"First three decay times (ms): {1e3/rates[:3]}")
    print(
        f"Eigenpair residual: {residual:.3g}; M-orthogonality error: {orthogonality:.3g}"
    )
    snapshots = np.array([0.02, 0.06, 0.18])
    times = np.unique(np.r_[np.linspace(0, 0.45, 241), snapshots, PULSE_DURATION])
    amplitudes = modal_response(
        rates,
        vectors.T @ coupling,
        current_derivative,
        times,
        max_step=PULSE_DURATION / 40,
    )
    coefficients = vectors @ amplitudes
    # Integrate the original coefficient ODE independently, without eigenvectors.
    factor = cho_factor(magnetic)
    evolution = -cho_solve(factor, resistance)
    forcing = -cho_solve(factor, coupling)
    direct = solve_ivp(
        lambda t, x: evolution @ x + forcing * current_derivative(t),
        (times[0], times[-1]),
        np.zeros(basis.n_dofs),
        method="BDF",
        jac=evolution,
        t_eval=times,
        rtol=1e-7,
        atol=1e-11,
        max_step=PULSE_DURATION / 40,
    )
    if not direct.success:
        raise RuntimeError(direct.message)
    difference = coefficients - direct.y
    error = np.sqrt(
        np.max(np.sum(difference * (magnetic @ difference), axis=0))
        / np.max(np.sum(coefficients * (magnetic @ coefficients), axis=0))
    )
    print(f"Modal versus direct transient relative error in M norm: {error:.3g}")
    if error > 1e-5:
        raise RuntimeError("modal/direct transient mismatch")
    stop = np.searchsorted(times, PULSE_DURATION)
    free = amplitudes[:, stop, None] * np.exp(
        -rates[:, None] * (times[None, stop:] - PULSE_DURATION)
    )
    decay_error = np.linalg.norm(amplitudes[:, stop:] - free) / np.linalg.norm(free)
    energy = np.sum(amplitudes**2, axis=0) / 2
    if decay_error > 1e-5 or np.any(np.diff(energy[stop:]) > 1e-10 * energy.max()):
        raise RuntimeError("free-decay verification failed")
    fields = reconstruct_currents(
        basis, coefficients[:, np.searchsorted(times, snapshots)]
    )
    # Match physical outward fluxes across each face, including the exterior.
    faces = mesh.vertices[mesh.tetrahedra[:, np.asarray(TETRA_FACES)]]
    areas = (
        np.cross(faces[:, :, 1] - faces[:, :, 0], faces[:, :, 2] - faces[:, :, 0]) / 2
    )
    flux = np.einsum("cfi,cit->cft", areas, fields)
    balances = np.zeros((len(topology.faces), len(snapshots)))
    np.add.at(balances, topology.cell_faces.ravel(), flux.reshape(-1, len(snapshots)))
    leakage = np.max(np.abs(balances)) / np.max(np.abs(flux))
    if leakage > 1e-10:
        raise RuntimeError("current continuity/boundary flux check failed")
    print(
        f"Free-decay relative error: {decay_error:.3g}; normalized face-flux imbalance: {leakage:.3g}"
    )
    print(f"Peak stored magnetic energy: {energy.max():.6g} J")
    if args.check_convergence:
        refined_magnetic = inductance_matrix(
            basis, quadrature_order=3, touching_order=6
        )
        refined_rates = eigh(
            resistance, refined_magnetic, eigvals_only=True, subset_by_index=(0, 2)
        )
        polygon_refined = loop_coupling(basis, source_loop(192), quadrature_order=6)
        print(
            f"First-rate relative quadrature changes: {np.abs(refined_rates-rates[:3])/refined_rates}"
        )
        print(
            f"Loop polygon refinement change in modal coupling: {np.linalg.norm(vectors.T @ (polygon_refined-refined_coupling))/np.linalg.norm(vectors.T @ polygon_refined):.3g}"
        )
    plotter = make_plotter(
        basis,
        loop,
        rates,
        times,
        amplitudes,
        fields,
        snapshots,
        off_screen=bool(args.screenshot),
    )
    plotter.show(screenshot=args.screenshot)


if __name__ == "__main__":
    main()
