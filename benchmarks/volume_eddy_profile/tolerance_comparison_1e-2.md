# Slow modes at rtol=1e-2

The 808-cell, 649-coefficient sphere from `examples/volume/eddy_currents.py`
was assembled at `rtol=1e-2` and compared with a fresh `rtol=1e-4` reference.
Minimum order two, maximum order twelve, block size sixteen, spatial ordering,
mesh, basis, and resistivity were held fixed. Production defaults were not changed.

## Time and quadrature work

| Measurement | rtol=1e-4 | rtol=1e-3, preceding run | rtol=1e-2 |
| --- | ---: | ---: | ---: |
| Uninstrumented assembly | 37.127 s | 31.074 s | 25.060 s |
| Adaptive unordered pairs | 317,373 (97.10%) | 259,745 (79.47%) | 178,579 (54.64%) |
| Direct order-two unordered pairs | 0 | 0 | 1,068 |
| Direct order-three unordered pairs | 0 | 5,087 | 75,782 |
| Direct order-four unordered pairs | 9,463 | 62,004 | 71,407 |
| Directed pair evaluations at order six | 64,834 | 12,080 | 9 |
| Directed pair evaluations at order eight | 101 | 3 | 0 |
| Directed pair evaluations at order ten | 2 | 0 | 0 |

The fresh reference and loose assemblies are one sequential uninstrumented
observation each, giving **1.482x speedup / 32.5% less elapsed time**. The
preceding 1e-3 run used a 36.854 s reference, so its timing is provided as context
rather than as part of a repeated timing distribution.

A separate instrumented 1e-2 assembly reproduced the uninstrumented matrix at
`rtol=2e-11`, `atol=1e-20`. It evaluated orders two and four for 359,046 directed
pair instances; only nine (0.00251%) required order six, and all converged there.
These counts include reciprocity and repeated work in mixed diagonal blocks.
There are 356,350 unique directed adaptive pairs, nine of which end at order six.

In that instrumented run, orders two/four cost approximately 3.860/15.237 CPU
seconds; order six cost 0.00154 s. Near integration totals 20.699 s and direct
far integration 4.918 s. Inclusive CPU counters include instrumentation overhead;
use the separate uninstrumented runs for elapsed-time comparisons. Refinement
above order four is now negligible, while the initial order-four evaluation
remains the largest cost.

The new direct-rule q cutoffs are 0.235806 (order two), 0.364123 (order three),
and 0.455964 (order four), where `q=(r_T+r_U)/centroid_distance`. The adaptive
cutoff is therefore approximately `d < 2.193*(r_T+r_U)`.

## Slow modes and smooth current patterns

Decay rates and eigenvectors solve `R v = lambda M v` for resistivity 1e-8 ohm m.
All changes below are relative to the 1e-4 calculation on the same mesh, rather
than an independently established continuum solution.

Patterns are compared in the physical volume-current L2 inner product:
`R = rho * integral(J_i dot J_j) dV`. Multiplication by the transpose of a
Cholesky factor of R maps coefficient vectors into that metric. Principal angles
compare the resulting slow-mode subspaces and allow rotations within nearly
degenerate groups. Individual correlations additionally remove sign ambiguity.

| Slowest modes included | Maximum relative decay-rate change | Largest current-subspace angle | Minimum individual current correlation |
| --- | ---: | ---: | ---: |
| 3 | 1.21e-7 | 0.0000353 degrees | 0.999999999923 |
| 10 | 2.26e-7 | 0.001185 degrees | 0.999999999784 |
| 20 | 3.73e-7 | 0.001546 degrees | 0.999999986413 |

The maximum relative energy change restricted to the first twenty reference
modes is **5.47e-7**. Their largest principal-angle sine is **2.70e-5**, measuring
the largest relative component outside the corresponding loose-mode subspace.

## Global matrix changes

- Relative Frobenius matrix difference: **5.6613e-5**.
- Maximum absolute matrix-entry difference: **1.1971e-12 H**.
- Energy ratios across all coefficient vectors: **0.9994852754 to 1.0007335386**.
  The largest relative energy change is **0.0734%**.
- Maximum relative change among all 649 decay rates: **5.7272e-4**.
- Both matrices remain positive definite in the generalized eigensolves.

The global matrix and energy differences increase relative to the 1e-3 case,
even though the measured slow-rate differences happen to be smaller. Individual
modal differences need not change monotonically with quadrature tolerance.

For this mesh and the stated interest in slow, smooth modes, the results support
using 1e-2: slow rates and current patterns remain very close while assembly
time falls by about one third relative to 1e-4.

## Reproduction

```sh
.venv/bin/python benchmarks/compare_volume_eddy_tolerances.py --rtol 1e-2 --output benchmarks/volume_eddy_profile/tolerance_comparison_1e-2.json
```

The JSON report contains the raw timings, complete pair/refinement counts,
thresholds, first ten decay rates, energy ratios, and pattern comparisons.
