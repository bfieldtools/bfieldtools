# Cached integral dataclass

`bfieldtools.integrals_cached.TriangleIntegrals` is derived from the current
functional `bfieldtools.integrals` implementation. The functional API remains
independent. The new class holds R and tn and exposes lazy cached properties:

```python
from bfieldtools.integrals_cached import TriangleIntegrals

triangle = TriangleIntegrals(R, tn, edges=edges)  # edges are optional
potential = triangle.triangle_potential_uniform
solid_angle = triangle.omega
```

Edges and optional supplied triangle areas seed their respective caches. Without
supplied values they are derived lazily. Vertex distances are shared by gamma0
and omega; volume face heights reuse d_distance. All implemented triangle and
dipole formulas from the functional module are available as cached properties.
The existing x_distance2 placeholder has no implemented formula to migrate.

Input arrays and cached values must not be mutated while the instance is in use.
The dataclass prevents field reassignment, borrows supplied arrays without
copying, and does not change their writeability. Volume assembly creates and
releases one instance per evaluation chunk so caches cannot grow across chunks.

Five warm repetitions, with alternating execution order, compared the exact
functional implementation immediately before this refactor against the new
cached implementation:

| Workload | Functional | Cached | Speedup |
| --- | ---: | ---: | ---: |
| 256 points / 64 triangle faces | 1.58 ms | 1.49 ms | 1.06x |
| 1,024 near points / 162 sources | 82.96 ms | 76.83 ms | 1.08x |
| 48-cell / 25-DOF inductance | 329.83 ms | 303.85 ms | 1.09x |
| 162-cell / 109-DOF inductance | 3.745 s | 3.461 s | 1.08x |

All compared numerical outputs were identical. The default 162-cell assembly
peaked at 5.86 MiB of traced allocations; the previous functional measurement
was about 5.73 MiB. Process RSS includes the imported numerical/visualization
stack, allocator retention, and prepared mesh/basis and is reported separately.

Validation: 125 numerical tests passed, including all volume suites and new
comparisons of cached exact, planar, normalized, approximate, dipole, and vertex
formulas with the functional module. Additional checks cover supplied read-only
geometry, noncontiguous inputs, single triangles, custom regularization,
frozen inputs, lazy planar evaluation, and one shared vertex-distance evaluation.

Artifacts:

- `volume_cached_integrals_results.json`: timing samples and numerical errors.
- `volume_cached_integrals_memory_results.json`: fresh-process memory profiles.
- `volume_cached_line_profile.txt` and `.json`: updated cached-property profiling.

The benchmark runner loads functional helpers or the cached class from the same
baseline revision as the volume implementation. Saved cached volume baseline
files additionally accept `--cached-baseline-file`. Line profiling uses the
cached-property getter functions directly.
