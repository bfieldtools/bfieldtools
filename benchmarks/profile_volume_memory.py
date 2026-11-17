"""Measure volume integral memory in a fresh process per configuration.

Run from the repository root:
    .venv/bin/python benchmarks/profile_volume_memory.py --output benchmarks/volume_memory_results.json

NumPy allocations are included in tracemalloc's peak. RSS is sampled every 1 ms;
its baseline includes imports and mesh/basis preparation. Timings are collected
separately from memory profiling, using three warmed calls. No package defaults
are changed when profiling alternative chunk sizes.
"""

from __future__ import annotations

import argparse
from functools import partial
import json
from pathlib import Path
import platform
import subprocess
import sys
import threading
from time import perf_counter
import tracemalloc
from unittest.mock import patch

import numpy as np
import psutil

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bfieldtools.volume import integrals, impedance
from tests.test_volume_inductance import basis_for
from tests.test_volume_topology import cube_mesh


def worker(config):
    mesh = cube_mesh(config.get("mesh_n", 3))
    points = np.random.default_rng(20261117).random((1024, 3))
    if config["kind"] == "far":
        points += 100
    elif config["kind"] == "mixed":
        points[512:] += 100
    kernel = partial(
        integrals.tetrahedron_potential,
        point_chunk_size=config["point_chunk_size"],
        source_chunk_size=config["source_chunk_size"],
    )
    basis = basis_for(mesh) if config["kind"] == "inductance" else None
    if basis is not None:
        # Materialize cached sparse basis operators before the measurement.
        _ = basis.operators

        def calculate():
            with patch.object(integrals, "tetrahedron_potential", kernel):
                return impedance.inductance_matrix(
                    basis, block_size=config["block_size"]
                )

    else:

        def calculate():
            return kernel(mesh, points)

    # Initialize quadrature caches before sampling, but keep array allocators
    # cold by measuring memory before warming the calculation itself.
    for order in (4, 6, 8, 10, 12):
        integrals.tetrahedron_quadrature(order)
    process = psutil.Process()
    rss_before = process.memory_info().rss
    rss_peak = [rss_before]
    stop = threading.Event()

    def sample():
        while not stop.wait(0.001):
            rss_peak[0] = max(rss_peak[0], process.memory_info().rss)

    thread = threading.Thread(target=sample, daemon=True)
    thread.start()
    tracemalloc.start()
    output = calculate()
    live, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    rss_peak[0] = max(rss_peak[0], process.memory_info().rss)
    stop.set()
    thread.join()
    row = dict(config)
    row.update(
        ncells=len(mesh.tetrahedra),
        output_shape=list(output.shape),
        output_bytes=output.nbytes,
        traced_peak_bytes=peak,
        traced_live_after_call_bytes=live,
        rss_baseline_bytes=rss_before,
        rss_peak_bytes=rss_peak[0],
        rss_increase_bytes=rss_peak[0] - rss_before,
    )
    if basis is not None and (
        config["block_size"],
        config["point_chunk_size"],
        config["source_chunk_size"],
    ) != (16, 256, 16):
        reference = impedance.inductance_matrix(basis)
        np.testing.assert_allclose(output, reference, rtol=2e-11, atol=1e-20)
        row["relative_difference_from_defaults"] = float(
            np.linalg.norm(output - reference) / np.linalg.norm(reference)
        )
        del reference
    del output
    samples = []
    for _ in range(3):
        start = perf_counter()
        output = calculate()
        samples.append(perf_counter() - start)
        del output
    row["median_seconds"] = float(np.median(samples))
    row["samples_seconds"] = samples
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker")
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(json.loads(args.worker))))
        return
    configs = [
        dict(kind=kind, point_chunk_size=256, source_chunk_size=16)
        for kind in ("near", "far", "mixed")
    ]
    configs.extend(
        [
            dict(
                kind="inductance",
                mesh_n=2,
                block_size=16,
                point_chunk_size=256,
                source_chunk_size=16,
            ),
            dict(
                kind="inductance",
                mesh_n=3,
                block_size=16,
                point_chunk_size=256,
                source_chunk_size=16,
            ),
            dict(
                kind="inductance",
                mesh_n=3,
                block_size=64,
                point_chunk_size=1024,
                source_chunk_size=64,
            ),
        ]
    )
    report = {
        "platform": platform.platform(),
        "numpy": np.__version__,
        "revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "method": "fresh process per case; tracemalloc includes NumPy; 1 ms RSS sampling; separate median of 3 warm timings",
        "cases": [],
    }
    for config in configs:
        print(f"Running {config}...", flush=True)
        result = subprocess.run(
            [sys.executable, __file__, "--worker", json.dumps(config)],
            text=True,
            capture_output=True,
            check=True,
        )
        row = json.loads(result.stdout)
        report["cases"].append(row)
        mib = 1024**2
        print(
            f"  traced peak {row['traced_peak_bytes']/mib:.2f} MiB; "
            f"RSS {row['rss_baseline_bytes']/mib:.1f} -> {row['rss_peak_bytes']/mib:.1f} MiB; "
            f"{row['median_seconds']:.3f} s",
            flush=True,
        )
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
