# 808-cell eddy_currents assembly profile

The exact mesh and basis setup from `examples/volume/eddy_currents.py` was used:
the bundled sphere scaled to radius 0.1 m, 808 tetrahedra, and 649 tree–cotree
current coefficients. Only `inductance_matrix(basis)` was measured. Mesh/basis
construction, resistance, eigensolving, loop coupling, and plotting are excluded.
Settings: block size 16, minimum order 2, `rtol=1e-4`, maximum order 12, and
spatial ordering enabled. This profiles the current working-tree integrator,
including the spatial/far-rule changes, rather than the unmodified base commit.

On macOS arm64, Python 3.11.17, NumPy 2.4.6, SciPy 1.17.1:

- Uninstrumented assembly: **37.356 s wall / 37.317 s process CPU** (one run).
- Separate cProfile assembly: **42.111 s**, using the profiler's default timer.
- Separate dispatch/refinement instrumentation: **38.815 s wall**.
- Profiled and instrumented matrices both matched the uninstrumented matrix at
  `rtol=2e-11`, `atol=1e-20`.

The timings are a single observation per mode, not a repeated timing distribution.

## Pair dispatch

Counts below refer to unique unordered pairs, including self interactions.

| Dispatch | Cell pairs | Fraction |
| --- | ---: | ---: |
| Adaptive/analytic integration | 317,373 | 97.105% |
| Direct far quadrature, order 4 | 9,463 | 2.895% |
| Direct far quadrature, orders 2 or 3 | 0 | 0% |
| Total | 326,836 | 100% |

Of the 317,373 adaptive pairs, **285,100 have disjoint bounding spheres**.
Only 32,273 pairs have overlapping spheres, including the 808 self pairs.
The radius-based Coulomb-series bound is conservative: many geometrically
separated pairs still take the adaptive/analytic path.

There are 1,326 upper-triangular blocks: 673 contain only adaptive pairs and
653 mix adaptive and direct-far pairs. **No block is entirely direct-far.**
Mixed blocks dispatch row subgroups, producing 11,121 near calls, 22,193
directed near calls, and 2,889 direct-far calls. There are 59,044 calls to the
potential evaluator and 92,221 analytic triangle batches.

## Refinement

Every adaptive interaction evaluates orders two and four. Higher orders are
evaluated only when the preceding comparison exceeds `rtol`. The table counts
executed directed cell-pair evaluations, including reciprocity and small amounts
of repeated work in mixed diagonal blocks; these are not unique unordered pairs.

| Order | Directed pair evaluations | Still active after comparison | Integration CPU time |
| --- | ---: | ---: | ---: |
| 2 | 634,400 | Not checked yet | 4.337 s |
| 4 | 634,400 | 64,834 | 24.002 s |
| 6 | 64,834 | 101 | 7.765 s |
| 8 | 101 | 2 | 0.037 s |
| 10 | 2 | 0 | 0.003 s |

Thus **10.22% of executed directed pairs require order six**. Counting each
directed pair once gives 10.18% requiring refinement beyond order four.
No pair reaches order twelve or fails convergence.

Orders two and four account for **28.339 s** of the instrumented integration
CPU time; orders six and above account for **7.805 s**. These times include
the instrumented point/source classification inside the potential evaluator,
so use them as a work split rather than uninstrumented timing predictions.

Of 59,734,656 outer-quadrature-point/source-cell evaluations inside the adaptive
path, 59,293,057 (99.26%) use analytic triangle potentials and 441,599 use the
existing order-four source rule beyond eight source radii. This internal source
fallback is distinct from the direct double-volume far-pair dispatcher.

## Hot functions

Inclusive times below are from the separate cProfile run. Rows are nested and
must not be added together.

| Function | Inclusive time | Fraction of profiled assembly |
| --- | ---: | ---: |
| `_near_interactions` | 40.478 s | 96.12% |
| `_evaluate_potential` | 37.970 s | 90.17% |
| `TriangleIntegrals.triangle_potential_uniform` | 28.150 s | 66.85% |
| `TriangleIntegrals.gamma0` | 15.572 s | 36.98% |
| `TriangleIntegrals.omega` | 7.051 s | 16.74% |
| `TriangleIntegrals.x_distance` | 4.455 s | 10.58% |
| `_SourceGeometry.from_mesh` | 1.196 s | 2.84% |
| `_regular_interactions` | 0.566 s | 1.34% |
| `_far_orders` | 0.051 s | 0.12% |

Self time identifies NumPy `c_einsum` (11.486 s), `_evaluate_potential` itself
(6.778 s), and `np.roll` (5.805 s self / 6.615 s inclusive, 737,768 calls) as
the largest contributors. Triangle determinants add 1.864 s self time.
The assembly spends about 0.872 s outside its interaction iterator, including
basis contractions and sparse indexing. Geometry integration dominates.

## Scaling interpretation

Linear near-field work requires a bounded number of local partners per cell,
bounded element aspect ratios and local size ratios, and bounded quadrature work
per pair. Under uniform refinement, a fixed separation-to-cell-size threshold
eventually supplies that locality. Direct evaluation of the remaining far pairs
is quadratic; spatial sorting alone does not provide hierarchical compression.

The present adaptive fallback is much broader than touching or nearest-neighbor
interactions. At the example's resolution it covers 97.1% of all pairs, and
99.26% of its point/source evaluations still use analytic potentials. Thus the
analytical work is effectively global at this resolution.

As a count-only check, the same dispatch criterion was applied to uniformly
refined unit cubes (no assembly timings were taken for this sequence):

| Cells | Adaptive unordered pairs | Fraction of all pairs |
| --- | ---: | ---: |
| 48 | 1,176 | 100% |
| 162 | 13,203 | 100% |
| 384 | 73,920 | 100% |
| 750 | 281,445 | 99.94% |
| 1,296 | 823,440 | 97.98% |
| 3,072 | 3,671,046 | 77.77% |

The near count grows nearly quadratically through the smaller cases. These
counts do not establish an asymptotic runtime exponent, but show that treating
the current analytic path as linear at these sizes is inaccurate.

The strongest optimization targets are reducing conservative fallback for
separated pairs with a validated tighter bound or a calibrated middle-distance
rule, batching the many small near/refinement groups, and reducing the repeated
triangle-array permutations and contractions. Optimizing direct far quadrature
alone has little headroom in this particular example.

## Reproduction and artifacts

```sh
.venv/bin/python benchmarks/profile_volume_eddy_assembly.py
```

- `results.json`: settings, pair counts, order counts, separate timings, scaling
  classification counts, and complete function profile.
- `functions.txt`: cumulative and self-time tables.
- `assembly.prof`: raw profile, readable with `python -m pstats`.

No production integration code was changed for this profiling task.
