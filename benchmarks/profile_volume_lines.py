"""Line-level process-CPU timings for current volume assembly hot functions.

Uses Python tracing, with no optional profiler dependency. Line timings include
called functions, so do not sum timings across different functions. Separate
untraced timing measures overhead. This measures NumPy calls as whole operations,
not individual native instructions or memory-bandwidth stalls.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bfieldtools import integrals_cached as triangle
from bfieldtools.volume import integrals as volume
from bfieldtools.volume.impedance import inductance_matrix
from tests.test_volume_inductance import basis_for
from tests.test_volume_topology import cube_mesh


class CpuLines:
    def __init__(self, functions):
        self.codes = {fn.__code__ for fn in functions}
        self.frames = {}
        self.lines = defaultdict(lambda: defaultdict(lambda: [0, 0.0]))

    def trace(self, frame, event, arg):
        if event == "call":
            if frame.f_code not in self.codes:
                return None
            self.frames[frame] = (None, process_time())
            return self.trace
        if event not in ("line", "return") or frame not in self.frames:
            return self.trace
        now = process_time()
        previous, start = self.frames[frame]
        if previous is not None:
            self.lines[frame.f_code][previous][1] += now - start
        if event == "return":
            del self.frames[frame]
            return None
        self.lines[frame.f_code][frame.f_lineno][0] += 1
        self.frames[frame] = (frame.f_lineno, process_time())
        return self.trace

    def results(self, repeats):
        result = []
        for code, lines in self.lines.items():
            source = Path(code.co_filename).read_text().splitlines()
            total = sum(values[1] for values in lines.values())
            rows = [
                dict(
                    line=line,
                    hits=hits,
                    cpu_seconds=seconds / repeats,
                    percent_of_function=100 * seconds / total,
                    source=source[line - 1].strip(),
                )
                for line, (hits, seconds) in sorted(lines.items())
            ]
            result.append(
                dict(
                    function=code.co_name,
                    file=code.co_filename,
                    inclusive_cpu_seconds=total / repeats,
                    lines=rows,
                )
            )
        return sorted(
            result, key=lambda row: row["inclusive_cpu_seconds"], reverse=True
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--mesh-n", type=int, default=3)
    parser.add_argument(
        "--output", type=Path, default=Path("benchmarks/volume_line_profile.json")
    )
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    basis = basis_for(cube_mesh(args.mesh_n))
    targets = [
        triangle.TriangleIntegrals.gamma0.func,
        triangle.TriangleIntegrals.omega.func,
        triangle.TriangleIntegrals.x_distance.func,
        triangle.TriangleIntegrals.d_distance.func,
        triangle.TriangleIntegrals.triangle_potential_uniform.func,
        triangle.TriangleIntegrals.vertex_distances.func,
        triangle.TriangleIntegrals.edge_dotprods1.func,
        triangle.TriangleIntegrals.edge_dotprods2.func,
        triangle.TriangleIntegrals.scaled_distances1.func,
        triangle.TriangleIntegrals.scaled_distances2.func,
        triangle.norm,
        triangle.determinant,
        triangle.cross,
        volume._evaluate_potential,
        volume._CellIntegrals._analytic,
        volume._CellIntegrals._regular,
    ]
    reference = inductance_matrix(basis)
    start_cpu, start_wall = process_time(), perf_counter()
    value = inductance_matrix(basis)
    untraced_cpu, untraced_wall = (
        process_time() - start_cpu,
        perf_counter() - start_wall,
    )
    np.testing.assert_allclose(value, reference, rtol=2e-11, atol=1e-20)
    profiler = CpuLines(targets)
    samples = []
    for repeat in range(args.repeats):
        start_cpu, start_wall = process_time(), perf_counter()
        sys.settrace(profiler.trace)
        try:
            value = inductance_matrix(basis)
        finally:
            sys.settrace(None)
        samples.append(
            dict(
                cpu_seconds=process_time() - start_cpu,
                wall_seconds=perf_counter() - start_wall,
            )
        )
        np.testing.assert_allclose(value, reference, rtol=2e-11, atol=1e-20)
        print(
            f"Traced assembly {repeat + 1}/{args.repeats}: {samples[-1]['cpu_seconds']:.3f} CPU s",
            flush=True,
        )
    functions = profiler.results(args.repeats)
    report = dict(
        platform=platform.platform(),
        numpy=np.__version__,
        ncells=len(basis.mesh.tetrahedra),
        ndofs=basis.n_dofs,
        repeats=args.repeats,
        untraced_cpu_seconds=untraced_cpu,
        untraced_wall_seconds=untraced_wall,
        traced_samples=samples,
        method="Python line trace; process CPU; inclusive per-line timings averaged across assemblies; NumPy kernels timed as calls",
        functions=functions,
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    text = [
        f"{len(basis.mesh.tetrahedra)}-cell / {basis.n_dofs}-DOF assembly; untraced CPU {untraced_cpu:.3f} s",
        "Line times are inclusive averages per assembly; hits count all traced runs. Do not sum across functions.",
    ]
    for fn in functions:
        text.append(
            f"\n{fn['function']}: {fn['inclusive_cpu_seconds']:.4f} CPU s / assembly"
        )
        text.append(" line    hits     CPU ms    % function   source")
        for row in fn["lines"]:
            text.append(
                f"{row['line']:5d} {row['hits']:7d} {1000*row['cpu_seconds']:10.3f} {row['percent_of_function']:10.2f}   {row['source']}"
            )
    args.output.with_suffix(".txt").write_text("\n".join(text) + "\n")
    print(f"Saved {args.output} and {args.output.with_suffix('.txt')}", flush=True)


if __name__ == "__main__":
    main()
