# Sphere volume refinement at rtol=1e-2

The bundled 808-cell sphere was remeshed with TetGen 0.8.4 using its original
faceted boundary and a smaller requested maximum unit-coordinate cell volume,
0.004. The resulting mesh has **1,804 cells**, exceeding the requested doubling.
It is a finer TetGen remesh rather than a nested subdivision of the coarse mesh.
Additional boundary triangles subdivide the existing planar boundary geometry.
The physical radius remains 0.1 m and resistivity 1e-8 ohm m.

The exact refined vertices/connectivity are archived in `mesh.npz`. Mesh
validation, positive orientation, face incidence, and quality checks pass.
Coarse and refined volumes agree to numerical precision.

## Timing comparison

Settings for both meshes: `rtol=1e-2`, minimum order two, maximum order twelve,
block size sixteen, spatial ordering. Each time is one sequential uninstrumented
assembly observation, on the same environment as the preceding tolerance tests.
Mesh generation, basis construction, classification, resistance assembly, and
eigensolving are excluded from inductance assembly time.

| Measurement | Coarse | Refined |
| --- | ---: | ---: |
| Tetrahedra | 808 | 1,804 |
| Vertices | 226 | 404 |
| Current coefficients | 649 | 1,603 |
| Inductance assembly wall time | **25.046 s** | **93.020 s** |
| Inductance assembly process CPU | 25.008 s | 92.905 s |
| Mesh/topology/basis setup | 0.0078 s | 0.0134 s |
| Solve for twenty slowest modes | 0.0209 s | 0.1240 s |
| Dense inductance matrix storage | 3.21 MiB | 19.60 MiB |
| Minimum mean-ratio quality | 0.2614 | 0.1527 |
| Median mean-ratio quality | 0.8790 | 0.8806 |
| Volume | 0.004047044680 m3 | 0.004047044680 m3 |

Cells grow **2.233x**, current coefficients **2.470x**, total unique cell pairs
**4.981x**, and assembly time **3.714x**. The two timing points correspond to an
effective exponent of 1.634 with cell count. This is a comparison of these two
meshes, not evidence for an asymptotic runtime law.

## Quadrature work

Pair counts are unique unordered pairs including the diagonal.

| Dispatch | Coarse pairs | Refined pairs |
| --- | ---: | ---: |
| Adaptive/analytic | 178,579 (54.64%) | 582,841 (35.80%) |
| Direct order two | 1,068 | 59,570 |
| Direct order three | 75,782 | 650,540 |
| Direct order four | 71,407 | 335,159 |
| Total | 326,836 | 1,628,110 |

The analytic partner set grows **3.264x**, less than the full pair count, because
the smaller cells make more pairs eligible for the direct far rules. It still
grows faster than cell count at this resolution: the current analytic path has
not yet reached a bounded-neighborhood regime.

A separate refined instrumented run took 97.175 s and matched the uninstrumented
matrix at `rtol=2e-11`, `atol=1e-20`. It recorded:

- 1,168,458 executed directed cell-pair instances at each of orders two and four.
- 121 executed directed pairs needing order six (0.0104%); all converged there.
  Of 1,163,878 unique directed adaptive pairs, 107 required order six.
- No evaluations at orders eight or higher.
- 84,155,112 outer-point/source-cell potential evaluations, all analytic.
- Approximately 68.07 CPU seconds in near integration and 25.98 seconds in direct
  far integration. The near path remains about 70% of instrumented wall time.
- Order-two integration: 12.98 CPU seconds; order-four integration: 49.60 seconds;
  order-six refinement: 0.0213 seconds.

These inclusive CPU counters include nested instrumentation overhead and must
not be summed across functions. The coarse instrumented run from the preceding
experiment cost approximately 20.70 seconds in near integration and 4.92 seconds
in direct far integration. Near cost grows roughly with its pair count; the
larger far set increasingly contributes to total assembly time.

There are 6,441 block evaluations: 373 entirely adaptive, six entirely direct
order three, and 6,062 mixed. Most far work continues to be dispatched as row
subgroups within mixed blocks.

## Slow-mode changes

Sorted generalized eigenvalues solve `R v = lambda M v`. The benchmark computes
the twenty slowest eigenpairs separately from the assembly timer. Both magnetic
matrices pass the positive-definiteness check implicit in the generalized
eigensolver. Relative eigenpair residuals are 4.29e-14 (coarse) and 1.60e-13
(refined); M-orthogonality errors are below 7.2e-15.

| Decay-rate statistic | Coarse | Refined | Relative change |
| --- | ---: | ---: | ---: |
| Slowest rate | 8.273180 s-1 | 8.167744 s-1 | -1.274% |
| Second rate | 8.279232 s-1 | 8.170674 s-1 | -1.311% |
| Third rate | 8.283092 s-1 | 8.171461 s-1 | -1.348% |
| Mean of slowest triplet | 8.278501 s-1 | 8.169960 s-1 | **-1.311%** |

The maximum relative change in the first ten sorted rates is 4.99%; in the first
twenty it is 7.54%. Coefficient eigenvectors from different meshes were not
compared directly because they represent different discrete current spaces.

The mesh-induced changes in this comparison are much larger than the quadrature
tolerance changes measured on the 808-cell mesh. Both runs use the same loose
quadrature setting, however, and no tighter-tolerance reference was assembled on
the refined mesh. Two meshes do not establish spatial convergence. This is
volume refinement of the same faceted domain, rather than refinement toward an
ideal spherical boundary.

## Reproduction

```sh
.venv/bin/python benchmarks/benchmark_volume_sphere_refinement.py
```

The script loads the archived refined mesh if present. Otherwise it remeshes the
original boundary, decreasing the volume target until at least twice the cells
are obtained, and saves the result. `results.json` contains the settings,
metadata, classifications, separate timings, all twenty slow rates, residuals,
and detailed refined-mesh instrumentation.
