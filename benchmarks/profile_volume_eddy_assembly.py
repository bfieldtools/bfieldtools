"""Profile the fixed-rule 808-cell eddy-current assembly, excluding setup/modes.

Run .venv/bin/python benchmarks/profile_volume_eddy_assembly.py.
Uninstrumented timing, cProfile, and path counters are separate assemblies.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from contextlib import ExitStack
import cProfile
import io
import json
from pathlib import Path
import platform
import pstats
import sys
from time import perf_counter, process_time
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bfieldtools.volume import (
    TetraMesh,
    TetraTopology,
    TreeCotreeGauge,
    build_edge_current_basis,
    inductance_matrix,
    load_example_tetra_mesh,
)
from bfieldtools.volume import integrals

SETTINGS = dict(block_size=16, quadrature_order=2, touching_order=4, spatial_order=True)


def classify_pairs(mesh):
    """Unique unordered contact counts; no global dense adjacency is retained."""
    counts = np.zeros(5, dtype=np.int64)
    columns = np.arange(len(mesh.tetrahedra))
    for start in range(0, len(columns), 64):
        rows = columns[start : start + 64]
        shared = np.sum(
            mesh.tetrahedra[rows, None, :, None]
            == mesh.tetrahedra[None, columns, None, :],
            axis=(2, 3),
        )
        counts += np.bincount(shared[rows[:, None] <= columns[None]], minlength=5)
    return dict(
        ncells=len(columns),
        unique_pairs=int(counts.sum()),
        direct_pairs=int(counts[0]),
        analytic_pairs=int(counts[1:].sum()),
        contacts=dict(
            zip(("non_touching", "vertex", "edge", "face", "self"), map(int, counts))
        ),
    )


def instrumented_assembly(basis, settings=None):
    """Separate kernel CPU/call counters; logical pair counts include block diagonals."""
    counts = defaultdict(lambda: dict(calls=0, cell_pairs=0, cpu_seconds=0.0))

    def timed(name):
        original = getattr(integrals._CellIntegrals, name)

        def call(geometry, rows, columns):
            start = process_time()
            value = original(geometry, rows, columns)
            row = counts[name]
            row["cpu_seconds"] += process_time() - start
            row["calls"] += 1
            row["cell_pairs"] += len(rows) * len(columns)
            return value

        return call

    with ExitStack() as stack:
        for name in ("_regular", "_analytic"):
            stack.enter_context(
                patch.object(integrals._CellIntegrals, name, timed(name))
            )
        start = perf_counter()
        result = inductance_matrix(
            basis, **(SETTINGS if settings is None else settings)
        )
        elapsed = perf_counter() - start
    return result, dict(wall_seconds=elapsed, functions=dict(counts))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmarks/volume_eddy_profile/fixed_profile"),
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    unit = load_example_tetra_mesh()
    mesh = TetraMesh(0.1 * unit.vertices, unit.tetrahedra)
    topology = TetraTopology.from_mesh(mesh)
    basis = build_edge_current_basis(
        mesh, topology, TreeCotreeGauge.from_topology(topology)
    )
    print("Running fixed-rule assembly...", flush=True)
    start = perf_counter()
    reference = inductance_matrix(basis, **SETTINGS)
    elapsed = perf_counter() - start
    report = dict(
        platform=platform.platform(),
        settings=SETTINGS,
        ncells=len(mesh.tetrahedra),
        ndofs=basis.n_dofs,
        pair_classification=classify_pairs(mesh),
        unprofiled_wall_seconds=elapsed,
    )
    profiler = cProfile.Profile()
    value = profiler.runcall(inductance_matrix, basis, **SETTINGS)
    np.testing.assert_allclose(value, reference, rtol=2e-11, atol=1e-20)
    profiler.dump_stats(str(args.output_dir / "assembly.prof"))
    output = io.StringIO()
    stats = pstats.Stats(profiler, stream=output)
    stats.strip_dirs().sort_stats("cumulative").print_stats(50)
    stats.sort_stats("tottime").print_stats(50)
    (args.output_dir / "functions.txt").write_text(output.getvalue().rstrip() + "\n")
    value, report["instrumentation"] = instrumented_assembly(basis)
    np.testing.assert_allclose(value, reference, rtol=2e-11, atol=1e-20)
    (args.output_dir / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
