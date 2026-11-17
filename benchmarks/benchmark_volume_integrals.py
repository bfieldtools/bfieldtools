"""Compare vectorized volume integrals with a previous implementation.

Run from the repository root, for example:
    .venv/bin/python benchmarks/benchmark_volume_integrals.py --baseline-ref HEAD

Use --baseline-file to compare against a saved integrals.py instead. Every case
is warmed up, checked numerically, and timed with alternating execution order.
"""

from __future__ import annotations

import argparse
from functools import partial
import json
import platform
from pathlib import Path
import subprocess
import sys
from time import perf_counter
import types
from unittest.mock import patch

import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bfieldtools import integrals as functional
from bfieldtools.integrals_cached import TriangleIntegrals
from bfieldtools.volume import integrals as current
from bfieldtools.volume import impedance
from bfieldtools.volume.mesh import TETRA_FACES
from tests.test_volume_inductance import basis_for
from tests.test_volume_topology import cube_mesh


def assemble_inductance(integrals, basis):
    # impedance imports the iterator directly. Replace that reference so the
    # complete assembly uses the selected implementation, including refinement.
    with patch.object(
        impedance,
        "iter_tetrahedron_interactions",
        integrals.iter_tetrahedron_interactions,
    ):
        return impedance.inductance_matrix(basis)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", default="HEAD")
    parser.add_argument("--baseline-file", type=Path)
    parser.add_argument("--surface-baseline-file", type=Path)
    parser.add_argument("--cached-baseline-file", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--only-inductance", action="store_true")
    parser.add_argument("--case", help="Run only cases containing this text")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.baseline_file:
        source = args.baseline_file.read_text()
        baseline_name = str(args.baseline_file)
    else:
        baseline_name = subprocess.check_output(
            ["git", "rev-parse", args.baseline_ref], text=True
        ).strip()
        source = subprocess.check_output(
            ["git", "show", f"{baseline_name}:bfieldtools/volume/integrals.py"],
            text=True,
        )
    baseline = types.ModuleType("bfieldtools.volume._benchmark_baseline")
    baseline.__package__ = "bfieldtools.volume"
    exec(compile(source, baseline_name, "exec"), baseline.__dict__)
    # Shared triangle helpers may also change. Load them from the same revision
    # so the baseline volume kernel cannot silently use the new implementation.
    if args.surface_baseline_file:
        surface_source = args.surface_baseline_file.read_text()
    elif not args.baseline_file:
        surface_source = subprocess.check_output(
            ["git", "show", f"{baseline_name}:bfieldtools/integrals.py"], text=True
        )
    else:
        surface_source = None
    if surface_source is not None:
        surface_baseline = types.ModuleType("bfieldtools._benchmark_surface_baseline")
        exec(compile(surface_source, baseline_name, "exec"), surface_baseline.__dict__)
        baseline.triangle_potential_uniform = (
            surface_baseline.triangle_potential_uniform
        )
    else:
        surface_baseline = functional
    if "TriangleIntegrals" in baseline.__dict__:
        if args.cached_baseline_file:
            cached_source = args.cached_baseline_file.read_text()
        elif not args.baseline_file:
            cached_source = subprocess.check_output(
                ["git", "show", f"{baseline_name}:bfieldtools/integrals_cached.py"],
                text=True,
            )
        else:
            parser.error(
                "A cached volume baseline file also requires --cached-baseline-file"
            )
        cached_baseline = types.ModuleType("bfieldtools._benchmark_cached_baseline")
        sys.modules[cached_baseline.__name__] = cached_baseline
        exec(compile(cached_source, baseline_name, "exec"), cached_baseline.__dict__)
        baseline.TriangleIntegrals = cached_baseline.TriangleIntegrals

    rng = np.random.default_rng(20261117)
    small, mesh = cube_mesh(1), cube_mesh(3)
    near = rng.random((256, 3))
    mixed = np.vstack((near[:128], near[128:] + [10, 0, 0]))
    faces = np.ascontiguousarray(
        mesh.vertices[mesh.tetrahedra[:16]][:, np.asarray(TETRA_FACES)]
    ).reshape(-1, 3, 3)
    triangle_normals = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0])
    triangle_normals /= np.linalg.norm(triangle_normals, axis=-1)[:, None]
    R = near[:, None, None] - faces[None]
    cases = [
        (
            "triangles: 256 points, 64 faces",
            "triangle_potential_uniform",
            (R, triangle_normals),
            {},
        ),
        (
            "triangles: 128 x 64 packed point/face pairs",
            "triangle_potential_uniform",
            (
                R[:128].reshape(1, -1, 3, 3),
                np.broadcast_to(
                    triangle_normals, (128, *triangle_normals.shape)
                ).reshape(-1, 3),
            ),
            {},
        ),
        (
            "potential: 1 point, 1 source",
            "tetrahedron_potential",
            (small, near[:1], [0]),
            {},
        ),
        (
            "potential: 1 point, 162 sources",
            "tetrahedron_potential",
            (mesh, near[:1]),
            {},
        ),
        (
            "potential: 256 near points, 162 sources",
            "tetrahedron_potential",
            (mesh, near),
            {},
        ),
        (
            "potential: 256 far points, 162 sources",
            "tetrahedron_potential",
            (mesh, near + 100),
            {},
        ),
        (
            "potential: 256 mixed points, 162 sources",
            "tetrahedron_potential",
            (mesh, mixed),
            {},
        ),
        (
            "potential: 1024 near points, 162 sources",
            "tetrahedron_potential",
            (mesh, rng.random((1024, 3))),
            {},
        ),
        (
            "interaction: 16 x 16 block",
            "tetrahedron_interaction_block",
            (mesh, np.arange(16), np.arange(16)),
            {},
        ),
        (
            "interaction: 16 x 16 off-diagonal block",
            "tetrahedron_interaction_block",
            (mesh, np.arange(16), np.arange(64, 80)),
            {},
        ),
        (
            "interaction: 48 x 48 matrix",
            "tetrahedron_interaction_matrix",
            (cube_mesh(2),),
            {},
        ),
    ]
    cases.extend(
        (
            f"inductance: {len(basis.mesh.tetrahedra)} cells, {basis.n_dofs} DOFs",
            "inductance_matrix",
            (basis,),
            {},
        )
        for basis in (basis_for(cube_mesh(2)), basis_for(mesh))
    )
    if args.only_inductance:
        cases = [case for case in cases if case[1] == "inductance_matrix"]
    if args.case:
        cases = [case for case in cases if args.case in case[0]]
        if not cases:
            parser.error("No benchmark cases match --case")
    report = {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "baseline": baseline_name,
        "repeats": args.repeats,
        "point_chunk_size": 256,
        "source_chunk_size": 16,
        "implementation": "bfieldtools.integrals_cached.TriangleIntegrals",
        "timing": "median wall time; warmup; alternating order",
        "inductance_settings": {
            "block_size": 16,
            "quadrature_order": 4,
            "rtol": 1e-4,
            "max_order": 12,
            "precomputed_interactions": False,
            "basis_construction_timed": False,
        },
        "cases": [],
    }
    for name, function, positional, keywords in cases:
        if function == "triangle_potential_uniform":
            if "TriangleIntegrals" in baseline.__dict__:
                old = lambda R, tn: baseline.TriangleIntegrals(
                    R, tn
                ).triangle_potential_uniform
            else:
                old = surface_baseline.triangle_potential_uniform
            new = lambda R, tn: TriangleIntegrals(R, tn).triangle_potential_uniform
        elif function == "inductance_matrix":
            old = partial(assemble_inductance, baseline)
            new = partial(assemble_inductance, current)
        else:
            old, new = getattr(baseline, function), getattr(current, function)
        print(f"Running {name}...", flush=True)
        expected = old(*positional, **keywords)
        actual = new(*positional, **keywords)
        np.testing.assert_allclose(
            actual,
            expected,
            rtol=2e-11,
            atol=1e-20 if function == "inductance_matrix" else 0,
        )
        timings = [[], []]
        for repeat in range(args.repeats):
            for index in [0, 1] if repeat % 2 == 0 else [1, 0]:
                start = perf_counter()
                (old, new)[index](*positional, **keywords)
                timings[index].append(perf_counter() - start)
        old_time, new_time = map(lambda times: float(np.median(times)), timings)
        row = {
            "name": name,
            "baseline_seconds": old_time,
            "vectorized_seconds": new_time,
            "speedup": old_time / new_time,
            "output_shape": list(actual.shape),
            "relative_frobenius_difference": float(
                np.linalg.norm(actual - expected) / np.linalg.norm(expected)
            ),
            "max_absolute_difference": float(np.max(abs(actual - expected))),
            "max_relative_difference": float(
                np.max(
                    abs(actual[expected != 0] - expected[expected != 0])
                    / abs(expected[expected != 0])
                )
            ),
            "baseline_samples_seconds": timings[0],
            "vectorized_samples_seconds": timings[1],
        }
        report["cases"].append(row)
        print(
            f"{name}: {old_time * 1000:.2f} -> {new_time * 1000:.2f} ms; "
            f"{old_time / new_time:.2f}x; "
            f"relative difference {row['max_relative_difference']:.2g}",
            flush=True,
        )
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
