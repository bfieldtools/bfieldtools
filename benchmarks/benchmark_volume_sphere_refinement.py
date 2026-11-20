"""Compare the 808-cell eddy sphere with at least twice as many volume cells.

The original faceted boundary, physical radius, material, and quadrature settings
are fixed. Meshing, topology/basis construction, pair classification, resistance,
and eigensolving are timed separately from uninstrumented inductance assembly.
A separate refined-mesh run counts adaptive work and checks repeatability.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time

import numpy as np
import scipy
from scipy.linalg import eigh
import trimesh

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmarks.profile_volume_eddy_assembly import (
    SETTINGS,
    classify_pairs,
    instrumented_assembly,
)
from bfieldtools.volume import (
    TetraMesh,
    TetraTopology,
    TreeCotreeGauge,
    build_edge_current_basis,
    inductance_matrix,
    load_example_tetra_mesh,
    resistance_matrix,
    tetrahedralize,
)


def refined_mesh(unit, folder, minimum_ratio):
    path = folder / "mesh.npz"
    if path.exists():
        with np.load(path, allow_pickle=False) as archive:
            fine = TetraMesh(archive["vertices"], archive["tetrahedra"])
            limit = float(archive["max_volume"])
        if len(fine.tetrahedra) < minimum_ratio * len(unit.tetrahedra):
            raise ValueError("Saved mesh has too few cells; use a new --output-dir")
        return fine, dict(
            method="loaded saved mesh; TetGen on original boundary",
            max_volume=limit,
            path=str(path),
        )
    boundary = trimesh.Trimesh(unit.vertices, unit.boundary_faces, process=False)
    boundary.remove_unreferenced_vertices()
    limit = 0.004
    attempts = []
    for _ in range(8):
        start = perf_counter()
        fine = tetrahedralize(boundary, max_volume=limit)
        attempts.append(
            dict(
                max_volume=limit,
                ncells=len(fine.tetrahedra),
                wall_seconds=perf_counter() - start,
            )
        )
        if len(fine.tetrahedra) >= minimum_ratio * len(unit.tetrahedra):
            np.savez_compressed(
                path,
                vertices=fine.vertices,
                tetrahedra=fine.tetrahedra,
                max_volume=limit,
            )
            return fine, dict(
                method="TetGen on original boundary",
                max_volume=limit,
                path=str(path),
                attempts=attempts,
            )
        limit *= 0.75
    raise RuntimeError("Could not reach the requested cell-count ratio")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmarks/volume_eddy_profile/fixed_rules_refinement"),
    )
    parser.add_argument("--minimum-cell-ratio", type=float, default=2.0)
    parser.add_argument("--quadrature-order", type=int, default=2)
    parser.add_argument("--touching-order", type=int, default=4)
    parser.add_argument(
        "--refined-mesh",
        type=Path,
        default=Path("benchmarks/volume_eddy_profile/refined_sphere/mesh.npz"),
    )
    args = parser.parse_args()
    if not np.isfinite(args.minimum_cell_ratio) or args.minimum_cell_ratio < 2:
        parser.error("--minimum-cell-ratio must be finite and at least two")
    if args.quadrature_order < 1 or args.touching_order < 1:
        parser.error("quadrature orders must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    unit = load_example_tetra_mesh()
    if args.refined_mesh.exists():
        with np.load(args.refined_mesh, allow_pickle=False) as archive:
            fine_unit = TetraMesh(archive["vertices"], archive["tetrahedra"])
        generation = dict(
            method="loaded archived refinement", path=str(args.refined_mesh)
        )
    else:
        fine_unit, generation = refined_mesh(
            unit, args.output_dir, args.minimum_cell_ratio
        )
    assert len(fine_unit.tetrahedra) >= args.minimum_cell_ratio * len(unit.tetrahedra)
    np.testing.assert_allclose(fine_unit.volumes.sum(), unit.volumes.sum(), rtol=1e-12)
    settings = dict(
        SETTINGS,
        quadrature_order=args.quadrature_order,
        touching_order=args.touching_order,
    )
    report = dict(
        platform=platform.platform(),
        python=platform.python_version(),
        numpy=np.__version__,
        scipy=scipy.__version__,
        settings=settings,
        radius_m=0.1,
        resistivity_ohm_m=1e-8,
        mesh_generation=generation,
        timing_method="one uninstrumented sequential assembly per mesh; separate setup, eigensolve, and instrumented refined assembly",
        cases={},
    )

    def save():
        (args.output_dir / "results.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )

    fine_reference = fine_basis = None
    for name, unscaled in (("coarse", unit), ("refined", fine_unit)):
        start = perf_counter()
        mesh = TetraMesh(0.1 * unscaled.vertices, unscaled.tetrahedra)
        topology = TetraTopology.from_mesh(mesh)
        basis = build_edge_current_basis(
            mesh, topology, TreeCotreeGauge.from_topology(topology)
        )
        case = dict(
            ncells=len(mesh.tetrahedra),
            nvertices=len(mesh.vertices),
            ndofs=basis.n_dofs,
            boundary_faces=len(mesh.boundary_faces),
            volume_m3=float(mesh.volumes.sum()),
            min_quality=float(mesh.cell_quality.min()),
            median_quality=float(np.median(mesh.cell_quality)),
            volume_min_m3=float(mesh.volumes.min()),
            volume_max_m3=float(mesh.volumes.max()),
            setup_wall_seconds=perf_counter() - start,
            dense_inductance_MiB=8 * basis.n_dofs**2 / 2**20,
        )
        case["pair_classification"] = classify_pairs(mesh)
        report["cases"][name] = case
        save()
        print(f"{name}: {case}", flush=True)
        print(f"Assembling {name} inductance...", flush=True)
        wall, cpu = perf_counter(), process_time()
        matrix = inductance_matrix(basis, **settings)
        case["assembly_wall_seconds"] = perf_counter() - wall
        case["assembly_cpu_seconds"] = process_time() - cpu
        save()
        print(f'{name} assembly: {case["assembly_wall_seconds"]:.3f}s', flush=True)

        start = perf_counter()
        resistance = resistance_matrix(basis, report["resistivity_ohm_m"]).toarray()
        case["resistance_wall_seconds"] = perf_counter() - start
        start = perf_counter()
        rates, vectors = eigh(resistance, matrix, subset_by_index=(0, 19))
        case["slow_modes_wall_seconds"] = perf_counter() - start
        case["slowest_20_decay_rates_per_s"] = rates.tolist()
        case["slowest_20_generalized_residual"] = float(
            np.linalg.norm(resistance @ vectors - (matrix @ vectors) * rates)
            / np.linalg.norm(resistance @ vectors)
        )
        case["slowest_20_M_orthogonality_error"] = float(
            np.linalg.norm(vectors.T @ matrix @ vectors - np.eye(20))
        )
        assert rates.min() > 0 and case["slowest_20_generalized_residual"] < 1e-9
        assert case["slowest_20_M_orthogonality_error"] < 1e-9
        save()
        if name == "refined":
            fine_reference, fine_basis = matrix, basis

    coarse, fine = report["cases"]["coarse"], report["cases"]["refined"]
    cell_ratio = fine["ncells"] / coarse["ncells"]
    time_ratio = fine["assembly_wall_seconds"] / coarse["assembly_wall_seconds"]
    rate_change = (
        np.asarray(fine["slowest_20_decay_rates_per_s"])
        / coarse["slowest_20_decay_rates_per_s"]
        - 1
    )
    report["comparison"] = dict(
        cell_ratio=cell_ratio,
        dof_ratio=fine["ndofs"] / coarse["ndofs"],
        pair_ratio=fine["pair_classification"]["unique_pairs"]
        / coarse["pair_classification"]["unique_pairs"],
        assembly_time_ratio=time_ratio,
        two_point_time_exponent=float(np.log(time_ratio) / np.log(cell_ratio)),
        relative_slowest_20_decay_rate_changes=rate_change.tolist(),
        max_relative_slowest_20_decay_rate_change=float(np.max(np.abs(rate_change))),
        mean_slowest_triplet_coarse=float(
            np.mean(coarse["slowest_20_decay_rates_per_s"][:3])
        ),
        mean_slowest_triplet_refined=float(
            np.mean(fine["slowest_20_decay_rates_per_s"][:3])
        ),
    )
    save()
    print(f'Comparison: {report["comparison"]}', flush=True)
    print("Instrumenting refined assembly...", flush=True)
    matrix, fine["instrumentation"] = instrumented_assembly(fine_basis, settings)
    np.testing.assert_allclose(matrix, fine_reference, rtol=2e-11, atol=1e-20)
    save()
    print(f"Saved {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
