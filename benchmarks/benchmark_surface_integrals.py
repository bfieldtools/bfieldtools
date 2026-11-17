"""Compare functional and cached integrals on surface ('2D') problems.

Run from the repository root:
    .venv/bin/python benchmarks/benchmark_surface_integrals.py --output benchmarks/surface_integrals_results.json

The production surface API is unchanged. Full self-inductance benchmarks replace
only the analytic triangle kernel within a context manager. Mesh preparation is
excluded; quadrature, near/far splitting, analytic diagonal correction, and
projection to vertex stream functions are included in assembly timings.
"""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from functools import partial
import io
import json
from pathlib import Path
import platform
import sys
from time import perf_counter
from unittest.mock import patch

import numpy as np
import scipy
import trimesh

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bfieldtools import integrals, mesh_impedance, mesh_magnetics
from bfieldtools.integrals_cached import TriangleIntegrals
from bfieldtools.quadratures import get_quad_points
from bfieldtools.utils import load_example_mesh


def cached_uniform(R, tn, planar=False):
    return TriangleIntegrals(R, tn, planar=planar).triangle_potential_uniform


def cached_dipole(R, tn, ta):
    return TriangleIntegrals(R, tn, ta=ta).triangle_potential_dipole_linear


def assemble(mesh, cached, **kwargs):
    kernel = cached_uniform if cached else integrals.triangle_potential_uniform
    with patch.object(mesh_magnetics, "triangle_potential_uniform", kernel):
        with redirect_stdout(io.StringIO()):
            return mesh_impedance.self_inductance_matrix(mesh, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--case", help="Run only cases containing this text")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    meshes = {
        "disc": load_example_mesh("unit_disc"),
        "plate": load_example_mesh("10x10_plane"),
        "sphere": trimesh.creation.icosphere(subdivisions=2, radius=1),
    }
    rng = np.random.default_rng(20261117)
    cases = []
    for label, mesh in meshes.items():
        planar = label != "sphere"
        triangles, normals, areas = (
            mesh.vertices[mesh.faces],
            mesh.face_normals,
            mesh.area_faces,
        )
        sampled = rng.integers(len(triangles), size=256)
        barycentric = rng.dirichlet(np.ones(3), size=256)
        points = np.einsum("pi,pij->pj", barycentric, triangles[sampled])
        R = points[:, None, None] - triangles[None]
        meta = dict(
            mesh=label,
            nfaces=len(mesh.faces),
            nvertices=len(mesh.vertices),
            npoints=256,
        )
        cases.append(
            (
                f"{label}: uniform kernel ({'planar' if planar else 'curved'})",
                partial(integrals.triangle_potential_uniform, R, normals, planar),
                partial(cached_uniform, R, normals, planar),
                dict(meta, kind="kernel", planar=planar),
            )
        )
        # Nonplanar evaluations use points offset from the source triangles.
        offset = points + 0.15 * np.linalg.norm(mesh.extents) * normals[sampled]
        R_offset = offset[:, None, None] - triangles[None]
        if planar:
            cases.append(
                (
                    f"{label}: uniform kernel (off-plane)",
                    partial(integrals.triangle_potential_uniform, R_offset, normals),
                    partial(cached_uniform, R_offset, normals),
                    dict(meta, kind="kernel", planar=False),
                )
            )
        cases.append(
            (
                f"{label}: linear dipole kernel",
                partial(
                    integrals.triangle_potential_dipole_linear, R_offset, normals, areas
                ),
                partial(cached_dipole, R_offset, normals, areas),
                dict(meta, kind="dipole"),
            )
        )
        weights, _ = get_quad_points(mesh.vertices, mesh.faces, ("dunavant", 3))
        nchunks = max(1, int(np.ceil(len(mesh.faces) * len(weights) / 128)))
        for approx_far in (False, True):
            settings = dict(
                Nchunks=nchunks,
                quad_degree=3,
                approx_far=approx_far,
                margin=2,
                chunk_clusters=False,
                planar=planar,
                analytic_self_coupling=True,
            )
            cases.append(
                (
                    f"{label}: self-inductance ({'approx far' if approx_far else 'exact'})",
                    partial(assemble, mesh, False, **settings),
                    partial(assemble, mesh, True, **settings),
                    dict(
                        meta,
                        kind="assembly",
                        npoints=len(mesh.faces) * len(weights),
                        **settings,
                    ),
                )
            )
    if args.case:
        cases = [case for case in cases if args.case in case[0]]
        if not cases:
            parser.error("No benchmark cases match --case")
    report = dict(
        platform=platform.platform(),
        python=platform.python_version(),
        numpy=np.__version__,
        scipy=scipy.__version__,
        repeats=args.repeats,
        method="Warm both implementations; median wall time; alternating execution order; new cache per evaluation call",
        surface_production_api_changed=False,
        cases=[],
    )
    for name, old, new, metadata in cases:
        print(f"Running {name}...", flush=True)
        expected, actual = old(), new()
        scale = float(np.max(abs(expected)))
        np.testing.assert_allclose(
            actual, expected, rtol=2e-11, atol=max(scale * 2e-14, 1e-25)
        )
        timings = [[], []]
        for repeat in range(args.repeats):
            for index in [0, 1] if repeat % 2 == 0 else [1, 0]:
                start = perf_counter()
                value = (old, new)[index]()
                timings[index].append(perf_counter() - start)
                del value
        old_seconds, new_seconds = [float(np.median(values)) for values in timings]
        row = dict(
            metadata,
            name=name,
            functional_seconds=old_seconds,
            cached_seconds=new_seconds,
            speedup=old_seconds / new_seconds,
            output_shape=list(actual.shape),
            relative_frobenius_difference=float(
                np.linalg.norm(actual - expected) / np.linalg.norm(expected)
            ),
            max_absolute_difference=float(np.max(abs(actual - expected))),
            functional_samples_seconds=timings[0],
            cached_samples_seconds=timings[1],
        )
        report["cases"].append(row)
        print(
            f"  {old_seconds*1000:.2f} -> {new_seconds*1000:.2f} ms; "
            f"{row['speedup']:.2f}x; relative norm difference "
            f"{row['relative_frobenius_difference']:.2g}",
            flush=True,
        )
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
