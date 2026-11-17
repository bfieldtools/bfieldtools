"""Profile streamed inductance assembly after mesh/basis preparation.

Example:
    .venv/bin/python benchmarks/profile_volume_cpu.py --output-dir /private/tmp/volume-cpu

Writes process-CPU cProfile data, readable cumulative/self-time tables, and JSON.
Unprofiled timings and order-by-order integration counters are separate runs.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import cProfile
from functools import partial
import io
import json
from pathlib import Path
import platform
import pstats
import subprocess
import sys
from time import perf_counter, process_time
from unittest.mock import patch

import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bfieldtools.volume import integrals, impedance
from tests.test_volume_inductance import basis_for
from tests.test_volume_topology import cube_mesh


def run_case(basis, config, output_dir):
    kernel = partial(
        integrals.tetrahedron_potential,
        point_chunk_size=config["point_chunk_size"],
        source_chunk_size=config["source_chunk_size"],
    )

    def calculate():
        with patch.object(integrals, "tetrahedron_potential", kernel):
            return impedance.inductance_matrix(basis, block_size=config["block_size"])

    reference = calculate()
    samples = []
    for _ in range(3):
        wall, cpu = perf_counter(), process_time()
        result = calculate()
        samples.append(
            {"wall_seconds": perf_counter() - wall, "cpu_seconds": process_time() - cpu}
        )
        np.testing.assert_allclose(result, reference, rtol=2e-11, atol=1e-20)
    print(
        f"{config['name']}: warm CPU {np.median([s['cpu_seconds'] for s in samples]):.3f}s",
        flush=True,
    )

    profiler = cProfile.Profile(timer=process_time)
    profiler.enable()
    result = calculate()
    profiler.disable()
    np.testing.assert_allclose(result, reference, rtol=2e-11, atol=1e-20)
    stem = output_dir / config["name"]
    profiler.dump_stats(str(stem.with_suffix(".prof")))
    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream)
    stats.strip_dirs().sort_stats("cumulative").print_stats(50)
    stats.sort_stats("tottime").print_stats(50)
    stem.with_suffix(".txt").write_text(stream.getvalue().rstrip() + "\n")

    entries = []
    # Keep original paths in JSON; the human-readable table uses short paths.
    raw_stats = pstats.Stats(profiler)
    for (filename, line, name), (
        primitive,
        calls,
        own,
        cumulative,
        _,
    ) in raw_stats.stats.items():
        entries.append(
            dict(
                file=filename,
                line=line,
                function=name,
                calls=calls,
                primitive_calls=primitive,
                self_cpu_seconds=own,
                cumulative_cpu_seconds=cumulative,
            )
        )
    orders = defaultdict(
        lambda: dict(
            calls=0,
            cell_pairs=0,
            point_source_pairs=0,
            cpu_seconds=0.0,
            single_row_calls=0,
        )
    )
    original_integrate = integrals._integrate_sources

    def count_integrate(mesh, rows, columns, order):
        cpu = process_time()
        value = original_integrate(mesh, rows, columns, order)
        record = orders[order]
        record["calls"] += 1
        record["cell_pairs"] += len(rows) * len(columns)
        record["point_source_pairs"] += len(rows) * len(columns) * order**3
        record["single_row_calls"] += int(len(rows) == 1)
        record["cpu_seconds"] += process_time() - cpu
        return value

    with patch.object(integrals, "_integrate_sources", count_integrate):
        result = calculate()
    np.testing.assert_allclose(result, reference, rtol=2e-11, atol=1e-20)
    report = dict(config)
    report.update(
        ncells=len(basis.mesh.tetrahedra),
        ndofs=basis.n_dofs,
        timing_samples=samples,
        median_cpu_seconds=float(np.median([s["cpu_seconds"] for s in samples])),
        median_wall_seconds=float(np.median([s["wall_seconds"] for s in samples])),
        profiled_cpu_seconds=raw_stats.total_tt,
        functions=sorted(
            entries, key=lambda row: row["cumulative_cpu_seconds"], reverse=True
        ),
        quadrature_orders=dict(sorted(orders.items())),
    )
    return report, reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("benchmarks/volume_cpu")
    )
    parser.add_argument("--mesh-n", type=int, default=3)
    parser.add_argument("--only-default", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    basis = basis_for(cube_mesh(args.mesh_n))
    configs = [
        dict(name="default", block_size=16, point_chunk_size=256, source_chunk_size=16),
        dict(
            name="larger_batches",
            block_size=64,
            point_chunk_size=1024,
            source_chunk_size=64,
        ),
    ]
    if args.only_default:
        configs = configs[:1]
    report = {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "method": "3 warm unprofiled timings; separate process_time cProfile run; separate order instrumentation run; mesh/basis excluded",
        "integration_settings": dict(quadrature_order=4, rtol=1e-4, max_order=12),
        "cases": [],
    }
    reference = None
    for config in configs:
        print(f"Profiling {config}...", flush=True)
        row, result = run_case(basis, config, args.output_dir)
        if reference is None:
            reference = result
        else:
            np.testing.assert_allclose(result, reference, rtol=2e-11, atol=1e-20)
            row["relative_difference_from_defaults"] = float(
                np.linalg.norm(result - reference) / np.linalg.norm(reference)
            )
        report["cases"].append(row)
        (args.output_dir / "results.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )
        print(f"Saved {config['name']} profile", flush=True)


if __name__ == "__main__":
    main()
