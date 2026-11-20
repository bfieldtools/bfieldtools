"""Compare fixed quadrature orders on the eddy-current sphere meshes.

Run .venv/bin/python benchmarks/compare_volume_eddy_orders.py.
The default (direct two, touching four) is compared with direct three/touching
six through the same implementation. This is an order-sensitivity check, not an
adaptive or guaranteed-accuracy matrix assembler. Production defaults are unchanged.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
from scipy.linalg import eigh, subspace_angles

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmarks.profile_volume_eddy_assembly import classify_pairs
from bfieldtools.volume import (
    TetraMesh,
    TetraTopology,
    TreeCotreeGauge,
    build_edge_current_basis,
    inductance_matrix,
    load_example_tetra_mesh,
    resistance_matrix,
)


def compare_modes(reference, candidate, resistance, reference_rates, reference_modes):
    # Since the reference is SPD, these eigenvalues test candidate positive
    # definiteness and bound magnetic-energy ratios for every coefficient vector.
    ratios = eigh(candidate, reference, eigvals_only=True)
    result = dict(
        relative_frobenius_matrix_difference=float(
            np.linalg.norm(candidate - reference) / np.linalg.norm(reference)
        ),
        max_absolute_matrix_difference_H=float(np.max(np.abs(candidate - reference))),
        energy_ratio_min=float(ratios.min()),
        energy_ratio_max=float(ratios.max()),
        positive_definite=bool(ratios.min() > 0),
    )
    if not result["positive_definite"]:
        return result
    rates, modes = eigh(resistance, candidate, subset_by_index=(0, 19))
    relative_rates = np.abs(rates - reference_rates) / reference_rates
    result["slowest_20_decay_rates_per_s"] = rates.tolist()
    result["max_relative_slowest_3_rate_change"] = float(relative_rates[:3].max())
    result["max_relative_slowest_10_rate_change"] = float(relative_rates[:10].max())
    result["max_relative_slowest_20_rate_change"] = float(relative_rates.max())
    result["generalized_residual"] = float(
        np.linalg.norm(resistance @ modes - (candidate @ modes) * rates)
        / np.linalg.norm(resistance @ modes)
    )
    assert rates.min() > 0 and result["generalized_residual"] < 1e-9
    factor = np.linalg.cholesky(resistance).T
    patterns = {}
    for count in (3, 10, 20):
        ref_fields = factor @ reference_modes[:, :count]
        fields = factor @ modes[:, :count]
        angles = subspace_angles(ref_fields, fields)
        ref_fields /= np.linalg.norm(ref_fields, axis=0)
        fields /= np.linalg.norm(fields, axis=0)
        projected = eigh(
            reference_modes[:, :count].T @ candidate @ reference_modes[:, :count],
            reference_modes[:, :count].T @ reference @ reference_modes[:, :count],
            eigvals_only=True,
        )
        patterns[str(count)] = dict(
            max_principal_angle_degrees=float(np.degrees(angles.max())),
            max_subspace_sine=float(np.sin(angles.max())),
            minimum_individual_current_correlation=float(
                np.min(np.abs(np.einsum("ij,ij->j", ref_fields, fields)))
            ),
            projected_energy_ratio_min=float(projected.min()),
            projected_energy_ratio_max=float(projected.max()),
            reference_relative_gap_to_next_mode=float(
                (reference_rates[count] - reference_rates[count - 1])
                / reference_rates[count - 1]
            )
            if count < 20
            else None,
        )
    result["current_patterns"] = patterns
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--meshes",
        choices=["coarse", "refined"],
        nargs="+",
        default=["coarse", "refined"],
    )
    parser.add_argument("--quadrature-order", type=int, default=2)
    parser.add_argument("--touching-order", type=int, default=4)
    parser.add_argument("--reference-quadrature-order", type=int, default=3)
    parser.add_argument("--reference-touching-order", type=int, default=6)
    parser.add_argument(
        "--refined-mesh",
        type=Path,
        default=Path("benchmarks/volume_eddy_profile/refined_sphere/mesh.npz"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmarks/volume_eddy_profile/fixed_order_comparison"),
    )
    args = parser.parse_args()
    if (
        min(
            args.quadrature_order,
            args.touching_order,
            args.reference_quadrature_order,
            args.reference_touching_order,
        )
        < 1
    ):
        parser.error("quadrature orders must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reference_settings = dict(
        quadrature_order=args.reference_quadrature_order,
        touching_order=args.reference_touching_order,
    )
    settings = dict(
        quadrature_order=args.quadrature_order, touching_order=args.touching_order
    )
    report = dict(
        platform=platform.platform(),
        reference_settings=reference_settings,
        settings=settings,
        timing_method="one assembly per mesh/order pair; setup and eigensolves excluded",
        cases={},
    )

    def save():
        (args.output_dir / "results.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )

    for name in args.meshes:
        if name == "coarse":
            unit = load_example_tetra_mesh()
        else:
            with np.load(args.refined_mesh, allow_pickle=False) as archive:
                unit = TetraMesh(archive["vertices"], archive["tetrahedra"])
        mesh = TetraMesh(0.1 * unit.vertices, unit.tetrahedra)
        topology = TetraTopology.from_mesh(mesh)
        basis = build_edge_current_basis(
            mesh, topology, TreeCotreeGauge.from_topology(topology)
        )
        resistance = resistance_matrix(basis, 1e-8).toarray()
        case = dict(
            ncells=len(mesh.tetrahedra), ndofs=basis.n_dofs, pairs=classify_pairs(mesh)
        )
        report["cases"][name] = case
        print(
            f"{name}: assembling reference fixed orders {reference_settings}",
            flush=True,
        )
        start = perf_counter()
        reference = inductance_matrix(basis, **reference_settings)
        case["reference_wall_seconds"] = perf_counter() - start
        rates, modes = eigh(resistance, reference, subset_by_index=(0, 19))
        print(f"{name}: assembling fixed orders {settings}", flush=True)
        start = perf_counter()
        matrix = inductance_matrix(basis, **settings)
        case["wall_seconds"] = perf_counter() - start
        np.testing.assert_array_equal(matrix, matrix.T)
        case["comparison"] = compare_modes(reference, matrix, resistance, rates, modes)
        save()
        assert case["comparison"]["positive_definite"]
        print(f"{name}: {case}", flush=True)
    print(f"Saved {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
