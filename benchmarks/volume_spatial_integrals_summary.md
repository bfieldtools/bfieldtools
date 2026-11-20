# Spatial volume integration

Baseline: `303f688eb3f9b65e5ef2c1a179568db094408844`.
Measurements use warmed calls, alternating baseline/current execution order,
and median wall time on macOS arm64, Python 3.11.17, NumPy 2.4.6, SciPy 1.17.1.
The baseline uses its original minimum order four; the current implementation
uses minimum order two. Both use `rtol=1e-4`, `max_order=12`, and block size 16.
Inductance timings share the current basis contraction and switch the complete
geometry integration implementation, including its blocking and order defaults.

| Case | Baseline | Current | Speedup | Relative Frobenius difference |
| --- | ---: | ---: | ---: | ---: |
| 16 × 16 local block | 34.04 ms | 31.45 ms | 1.08× | 6.27e-8 |
| 16 × 16 off-diagonal block | 66.12 ms | 53.87 ms | 1.23× | 1.26e-7 |
| 48 × 48 dense cell matrix | 304.80 ms | 289.70 ms | 1.05× | 1.05e-7 |
| 16 × 16 distant neighborhoods | 90.90 ms | 0.410 ms | 221.79× | 4.73e-11 |
| Inductance, 48 cells / 25 DOFs | 322.18 ms | 307.14 ms | 1.05× | 1.76e-7 |
| Inductance, 162 cells / 109 DOFs | 3657.88 ms | 2461.28 ms | 1.49× | 2.97e-7 |

The distant case takes two neighborhoods from a 162-cell unit-cube mesh,
separated by 30 coordinate units. Its speedup is specific to that separation;
it does not describe a full mesh assembly. The regular far rule evaluates
both volume integrals directly and needs neither triangle potentials nor
successive-order comparisons. Local interactions retain analytic inner
integration and empirical convergence checks.

Maximum entrywise relative cell-kernel differences were 9.05e-7 for the local
cases and 9.36e-11 for the distant case. For inductance, signed contractions
produce small entries by cancellation, so comparison uses the Frobenius norm.
The largest entrywise relative difference was 4.64e-4 in the 162-cell case;
the maximum absolute difference over that matrix was 1.22e-14 H.

Validation: 112 volume tests passed, covering independent integral references,
touching cells, positivity, energy, generalized modes, geometric scaling,
mixed near/far blocks, original-numbering scatter, permutation invariance,
far-order bounds against higher-order quadrature, and sliver rejection.

Reproduce from the repository root:

```sh
.venv/bin/python benchmarks/benchmark_volume_integrals.py --baseline-ref 303f688eb3f9b65e5ef2c1a179568db094408844 --case interaction --rtol 1e-4 --repeats 5 --output benchmarks/volume_spatial_integrals_results.json
.venv/bin/python benchmarks/benchmark_volume_integrals.py --baseline-ref 303f688eb3f9b65e5ef2c1a179568db094408844 --only-inductance --rtol 1e-4 --repeats 3 --output benchmarks/volume_spatial_inductance_results.json
```

Raw samples are in the corresponding JSON reports. The kernel bound and API
options are described in `bfieldtools/volume/INDUCTANCE_INTEGRALS.md`.
