"""Eddy currents in a solid torus driven by a coaxial current loop.

Run from the repository root::

    uv run --locked --extra meshing python examples/volume/torus_eddy_currents.py

The torus has major radius 0.1 m, tube radius 0.025 m, and resistivity 1e-8
ohm metres. TetGen fills a coarse faceted surface. A prescribed 1000 A smooth
40 ms pulse in a nearby coaxial loop excites both local and global currents.
All generalized modes are retained. After the pulse the currents decay freely.

The PyVista window shows the drive, net tube current, magnetic energy, and three
linked current snapshots. Colored cutaways show the actual cellwise current
magnitude; uniformly spaced vectors sample the current in their containing
tetrahedra. Arrow lengths are fixed to show direction; color carries magnitude.
``--show-mesh`` reveals the
cell edges. ``--no-show`` runs the numerical checks without opening a window;
``--screenshot torus.png`` saves the window off screen. Increase
``--major-sections``, ``--minor-sections``, and reduce ``--max-volume`` (m^3) for
mesh studies. ``--check-convergence`` tightens integration on the same mesh;
this does not establish spatial accuracy or skin-depth resolution. The
augmented basis uses a size-limited dense reference construction.
"""

import argparse

import numpy as np
import trimesh
from scipy.integrate import solve_ivp
from scipy.linalg import cho_factor, cho_solve, eigh

from bfieldtools.volume import (
    augmented_flux_basis,
    build_face_current_basis,
    build_topology,
    inductance_matrix,
    loop_coupling,
    modal_response,
    reconstruct_currents,
    resistance_matrix,
    tetrahedralize,
)
from bfieldtools.volume.mesh import TETRA_EDGES


PULSE_DURATION = 0.04
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


def source_loop(segments=96):
    angle = np.linspace(0, 2 * np.pi, segments + 1)
    vertices = np.column_stack(
        (0.15 * np.cos(angle), 0.15 * np.sin(angle), np.full_like(angle, 0.06))
    )
    vertices[-1] = vertices[0]
    return vertices


def cross_section_weights(mesh, angle):
    """Integrate constant cell currents across one radial half-plane exactly.

    An oblique plane avoids coincident mesh faces. Each intersected tetrahedron
    contributes its polygonal section area times the plane normal. This is an
    example-specific observable, not a general cut API.
    """
    radial = np.array([np.cos(angle), np.sin(angle), 0.0])
    normal = np.array([-np.sin(angle), np.cos(angle), 0.0])
    weights = np.zeros((len(mesh.tetrahedra), 3))
    tolerance = 1e-12 * np.ptp(mesh.vertices, axis=0).max()
    for cell, points in enumerate(mesh.vertices[mesh.tetrahedra]):
        distances = points @ normal
        if distances.min() >= 0 or distances.max() <= 0:
            continue
        if np.any(abs(distances) < tolerance):
            raise ValueError("cut passes through a mesh vertex; choose another angle")
        intersections = []
        for a, b in TETRA_EDGES:
            if distances[a] * distances[b] < 0:
                fraction = distances[a] / (distances[a] - distances[b])
                intersections.append(points[a] + fraction * (points[b] - points[a]))
        polygon = np.asarray(intersections)
        radii = polygon @ radial
        if radii.mean() < 0:
            continue
        if radii.min() <= 0:
            raise ValueError("the cut must intersect only one side of the torus")
        polygon = np.column_stack((radii, polygon[:, 2]))
        polygon -= polygon.mean(axis=0)
        polygon = polygon[np.argsort(np.arctan2(polygon[:, 1], polygon[:, 0]))]
        area = (
            abs(
                np.sum(
                    polygon[:, 0] * np.roll(polygon[:, 1], -1)
                    - polygon[:, 1] * np.roll(polygon[:, 0], -1)
                )
            )
            / 2
        )
        weights[cell] = area * normal
    return weights


def simulate(args):
    surface = trimesh.creation.torus(
        major_radius=0.1,
        minor_radius=0.025,
        major_sections=args.major_sections,
        minor_sections=args.minor_sections,
    )
    mesh = tetrahedralize(surface, max_volume=args.max_volume)
    topology = build_topology(mesh)
    space = augmented_flux_basis(topology)
    if space.n_global != 1:
        raise RuntimeError("a solid torus must have exactly one global current mode")
    basis = build_face_current_basis(mesh, space)
    print(
        f"{len(mesh.tetrahedra)} cells; {basis.n_dofs} currents = "
        f"{len(space.local_edges)} local + 1 global",
        flush=True,
    )
    resistance = resistance_matrix(basis, 1e-8).toarray()
    magnetic = inductance_matrix(
        basis,
        quadrature_order=args.quadrature_order,
        touching_order=args.touching_order,
    )
    rates, vectors = eigh(resistance, magnetic)
    residual = np.linalg.norm(
        resistance @ vectors - (magnetic @ vectors) * rates
    ) / np.linalg.norm(resistance @ vectors)
    orthogonality = np.linalg.norm(vectors.T @ magnetic @ vectors - np.eye(len(rates)))
    if rates[0] <= 0 or residual > 1e-10 or orthogonality > 1e-9:
        raise RuntimeError("generalized eigenpair verification failed")
    print(f"First three decay times (ms): {1e3 / rates[:3]}")
    loop = source_loop()
    coupling = loop_coupling(basis, loop)
    snapshots = np.array(
        [PULSE_DURATION / 4, 3 * PULSE_DURATION / 4, PULSE_DURATION + 1 / rates[0]]
    )
    times = np.unique(
        np.r_[
            np.linspace(0, PULSE_DURATION + 5 / rates[0], 241),
            snapshots,
            PULSE_DURATION,
        ]
    )
    amplitudes = modal_response(
        rates,
        vectors.T @ coupling,
        current_derivative,
        times,
        max_step=PULSE_DURATION / 40,
    )
    coefficients = vectors @ amplitudes

    # Independently integrate M dx/dt + R x = -m dI/dt in coefficient space.
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
    energy = np.sum(amplitudes**2, axis=0) / 2
    difference = coefficients - direct.y
    transient_error = np.sqrt(
        np.max(np.sum(difference * (magnetic @ difference), axis=0))
        / (2 * energy.max())
    )
    stop = np.searchsorted(times, PULSE_DURATION)
    free = amplitudes[:, stop, None] * np.exp(
        -rates[:, None] * (times[None, stop:] - PULSE_DURATION)
    )
    if (
        transient_error > 1e-5
        or np.linalg.norm(amplitudes[:, stop:] - free) > 1e-5 * np.linalg.norm(free)
        or np.any(np.diff(energy[stop:]) > 1e-10 * energy.max())
    ):
        raise RuntimeError("transient/free-decay verification failed")

    # Net current is independent of where the tube is cut. Local edge curls
    # contribute no net cut current; only the global coefficient contributes.
    cut_maps = []
    for angle in (0.137, 0.431):
        weights = cross_section_weights(mesh, angle)
        cut_maps.append(
            sum(
                operator.T @ weights[:, k] for k, operator in enumerate(basis.operators)
            )
        )
    np.testing.assert_allclose(cut_maps[0], cut_maps[1], atol=1e-11)
    np.testing.assert_allclose(cut_maps[0][:-1], 0, atol=1e-11)
    circulating = cut_maps[0] @ coefficients
    if np.max(abs(circulating)) < 1e-6:
        raise RuntimeError("the drive did not excite the global circulating current")
    fields = reconstruct_currents(
        basis, coefficients[:, np.searchsorted(times, snapshots)]
    )
    triangles = mesh.vertices[topology.faces]
    area_vectors = (
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        / 2
    )
    measured = np.einsum("cfi,cit->cft", area_vectors[topology.cell_faces], fields)
    fluxes = space.fluxes @ coefficients[:, np.searchsorted(times, snapshots)]
    leakage = np.max(abs(measured - fluxes[topology.cell_faces])) / np.max(
        abs(measured)
    )
    if leakage > 1e-10:
        raise RuntimeError(
            "physical face fluxes disagree with the conserved flux basis"
        )
    print(
        f"Peak circulating current: {np.max(abs(circulating)):.6g} A; "
        f"peak magnetic energy: {energy.max():.6g} J"
    )
    print(
        f"Modal/direct relative error: {transient_error:.3g}; "
        f"face-flux relative error: {leakage:.3g}"
    )
    if args.check_convergence:
        refined = inductance_matrix(
            basis,
            quadrature_order=args.quadrature_order + 1,
            touching_order=args.touching_order + 2,
        )
        refined_rates = eigh(
            resistance, refined, eigvals_only=True, subset_by_index=(0, 2)
        )
        refined_coupling = loop_coupling(basis, source_loop(192), quadrature_order=6)
        print(
            "First-rate changes with tighter quadrature:",
            abs(refined_rates - rates[:3]) / refined_rates,
        )
        print(
            "Modal coupling change with refined loop/integration:",
            np.linalg.norm(vectors.T @ (refined_coupling - coupling))
            / np.linalg.norm(vectors.T @ refined_coupling),
        )
    return basis, loop, times, circulating, energy, fields, snapshots


def make_plotter(
    basis,
    loop,
    times,
    circulating,
    energy,
    fields,
    snapshots,
    *,
    off_screen=False,
    show_mesh=False,
):
    """Linked magnitude cutaways and regularly sampled current directions."""
    import matplotlib.pyplot as plt
    import pyvista as pv

    plotter = pv.Plotter(shape=(2, 3), window_size=(1680, 1050), off_screen=off_screen)
    plotter.set_background("white")
    for panel, (values, title, ylabel) in enumerate(
        (
            (waveform(times), "Prescribed drive", "Loop current (A)"),
            (circulating, "Net tube current", "Circulating current (A)"),
            (energy, "Stored magnetic energy", "Energy (J)"),
        )
    ):
        figure, axes = plt.subplots(figsize=(5, 3.2))
        axes.plot(1e3 * times, values, color="#256b80", linewidth=2)
        axes.axvline(1e3 * PULSE_DURATION, color="0.6", linestyle=":")
        for t in snapshots:
            axes.axvline(1e3 * t, color="0.85", linewidth=0.8)
        axes.set(xlabel="Time (ms)", ylabel=ylabel, title=title)
        axes.grid(alpha=0.15)
        figure.tight_layout()
        plotter.subplot(0, panel)
        plotter.add_chart(pv.ChartMPL(figure))
        plt.close(figure)

    mesh = basis.mesh
    grid = mesh.to_pyvista()
    boundary = grid.extract_surface(algorithm="dataset_surface")
    maximum = np.linalg.norm(fields, axis=1).max()
    # Sample in physical space, independently of the tetrahedron sizes. Use the
    # same locations in every snapshot and retain the containing cell's vector.
    angles = np.linspace(0, 2 * np.pi, 24, endpoint=False) + 0.137
    samples = np.array(
        [
            [r * np.cos(angle), r * np.sin(angle), 0.008]
            for angle in angles
            for r in (0.09, 0.11)
        ]
    )
    cell_ids = grid.find_containing_cell(samples)
    samples, cell_ids = samples[cell_ids >= 0], cell_ids[cell_ids >= 0]
    for panel, t in enumerate(snapshots):
        plotter.subplot(1, panel)
        plotter.add_text(
            f"t = {1e3 * t:.1f} ms   |   drive = {float(waveform(t)):.0f} A",
            font_size=12,
            color="#183746",
        )
        plotter.add_text(
            "Color: current magnitude   |   Cyan arrows: direction",
            position="lower_left",
            font_size=9,
            color="#4a5960",
        )
        grid.cell_data["current_density"] = np.linalg.norm(fields[:, :, panel], axis=1)
        # Remove the upper front quarter: retain the complete back half and
        # the lower front half, exposing both surface and interior magnitudes.
        back = grid.clip(normal="y", origin=(0, 0, 0), invert=False)
        front = grid.clip(normal="y", origin=(0, 0, 0), invert=True)
        front = front.clip(normal="z", origin=(0, 0, 0), invert=True)
        cutaway = back.merge(front)
        plotter.add_mesh(boundary, color="#adc0c7", opacity=0.09)
        plotter.add_mesh(
            cutaway,
            scalars="current_density",
            preference="cell",
            cmap="inferno",
            clim=(0, maximum),
            lighting=False,
            show_edges=show_mesh,
            edge_color="#6b7380",
            edge_opacity=0.25,
            scalar_bar_args={
                "title": "|J| (A/m^2), shared scale",
                "position_x": 0.14,
                "position_y": 0.1,
                "width": 0.72,
                "height": 0.07,
                "title_font_size": 11,
                "label_font_size": 10,
                "fmt": "%.1e",
            },
        )
        vectors = fields[cell_ids, :, panel]
        lengths = np.linalg.norm(vectors, axis=1)
        active = lengths > maximum * 1e-10
        if active.any():
            points = pv.PolyData(samples[active])
            points["direction"] = vectors[active] / lengths[active, None]
            arrows = points.glyph(
                orient="direction",
                scale=False,
                factor=0.014,
                geom=pv.Arrow(tip_length=0.32, tip_radius=0.16, shaft_radius=0.05),
            )
            plotter.add_mesh(
                arrows,
                color="#4bd9e0",
                ambient=0.4,
                diffuse=0.6,
            )
        plotter.add_mesh(
            pv.lines_from_points(loop), color="#b54f4f", line_width=2, opacity=0.6
        )
        plotter.add_axes()
        plotter.view_vector((2, -3, 2.5))
        plotter.camera.zoom(1.15)
    plotter.link_views(views=(3, 4, 5))
    return plotter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--major-sections", type=int, default=12)
    parser.add_argument("--minor-sections", type=int, default=8)
    parser.add_argument("--max-volume", type=float)
    parser.add_argument("--quadrature-order", type=int, default=2)
    parser.add_argument("--touching-order", type=int, default=4)
    parser.add_argument("--check-convergence", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument(
        "--show-mesh",
        action="store_true",
        help="show tetrahedral edges on the cutaways",
    )
    parser.add_argument("--screenshot")
    args = parser.parse_args()
    if min(args.major_sections, args.minor_sections) < 6:
        parser.error("surface section counts must be at least six")
    result = simulate(args)
    if not args.no_show or args.screenshot:
        if args.screenshot:
            import matplotlib

            matplotlib.use("Agg")
        plotter = make_plotter(
            *result, off_screen=bool(args.screenshot), show_mesh=args.show_mesh
        )
        plotter.show(screenshot=args.screenshot)


if __name__ == "__main__":
    main()
