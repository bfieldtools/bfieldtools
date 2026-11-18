# Eddy-current assembly: rtol 1e-4 versus 1e-3

The same 808-cell sphere and 649-coefficient basis as `eddy_currents.py` were
assembled at both tolerances. All other integration settings were unchanged:
minimum order two, maximum order twelve, block size sixteen, spatial ordering.
Production defaults were not changed.

## Separation thresholds

The thresholds solve `2*(1+q)*q**(2*n)/(1-q) = rtol`, where
`q = (r_T+r_U)/centroid_distance`. The lowest eligible direct rule is chosen;
otherwise integration uses the adaptive/analytic path.

| Direct order | q cutoff at 1e-4 | q cutoff at 1e-3 |
| --- | ---: | ---: |
| 2 | 0.080755 | 0.139403 |
| 3 | 0.180604 | 0.257998 |
| 4 | 0.270550 | 0.352666 |

Thus the approximate distance cutoff for adaptive integration moves from
`3.696*(r_T+r_U)` to `2.836*(r_T+r_U)`.

## Observed work and time

| Measurement | rtol=1e-4 | rtol=1e-3 |
| --- | ---: | ---: |
| Adaptive unordered pairs | 317,373 (97.10%) | 259,745 (79.47%) |
| Direct order-two unordered pairs | 0 | 0 |
| Direct order-three unordered pairs | 0 | 5,087 |
| Direct order-four unordered pairs | 9,463 | 62,004 |
| Uninstrumented wall time | 36.854 s | 31.074 s |
| Executed directed near pairs requiring order six | 64,834 (10.22%) | 12,080 (2.32%) |
| Directed pair evaluations at order eight | 101 | 3 |
| Directed pair evaluations at order ten | 2 | 0 |

The timing observations are one sequential uninstrumented run per tolerance,
with the stricter run first: **1.186x speedup / 15.7% less elapsed time**.
They are not repeated timing distributions. Refinement counts at 1e-4 come from
the preceding detailed profile; counts at 1e-3 come from a new separate
instrumented run, which matched its uninstrumented matrix at `rtol=2e-11`,
`atol=1e-20`.

Refinement above order four costs 7.805 s of instrumented integration CPU at
1e-4 and 1.673 s at 1e-3. Direct far quadrature increases from approximately
0.550 s to 3.401 s. More mixed blocks also fragment the near evaluations:
directed near calls increase from 22,193 to 32,828 despite fewer near pairs.
These inclusive CPU counters contain instrumentation overhead and are useful
for explaining the work split, rather than as direct uninstrumented predictions.

## Changes relative to the 1e-4 result

- Relative Frobenius inductance-matrix difference: **1.6935e-5**.
- Maximum absolute matrix-entry difference: **3.5188e-13 H**.
- Magnetic energy ratio over all nonzero coefficient vectors lies between
  **0.9997518424 and 1.0002518177**, computed from the extremal generalized
  eigenvalues of `M_loose x = lambda M_reference x`. The largest relative
  energy change is therefore approximately **0.0252%**.
- Maximum relative change among the ten slowest decay rates: **6.5202e-7**.
- Maximum relative change among all decay rates: **2.1273e-4**.
- Both matrices remain positive definite, as checked by generalized symmetric
  eigensolvers.

These are differences from the 1e-4 calculation, not independently established
errors against the exact continuum problem. The mesh discretization is unchanged.
The looser tolerance gives a modest speedup with very small changes to the
slowest modes in this example; the analytic path still handles most pairs.

## Reproduction

```sh
.venv/bin/python benchmarks/compare_volume_eddy_tolerances.py
```

Full settings, classifications, timings, energy ratios, decay rates, and
instrumentation counters are in `tolerance_comparison.json`.
