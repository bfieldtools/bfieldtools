# Analytic self/touching pairs and fixed quadrature elsewhere

This experiment replaces the current distance/error-bound dispatch with one
topological rule: cells sharing any vertex use the existing analytic inner
integration and adaptive outer rule; all other pairs use fixed double-product
quadrature. Contact includes self, face, edge, and vertex neighbors. Two fixed
orders were tested: two (eight points per cell, 64 point pairs) and three
(27 points per cell, 729 point pairs).

The experiment is isolated in `benchmark_volume_touching_policy.py`, using a
scoped replacement of the private order selector. Production behavior and
defaults are unchanged. All implementations use the same integration kernels,
spatial blocks, streamed basis contractions, and mesh/basis arrays. The
experimental non-touching rules have **no a priori error bound or refinement**;
`rtol=1e-2` controls only analytic self/touching interactions. Other settings are
minimum order two, maximum order twelve, block size sixteen, and spatial ordering.

## Assembly timings

Each is one sequential assembly observation; setup, classification, and
comparison eigensolves are outside the timer. The refined sphere is the archived
1,804-cell mesh with the same faceted boundary and radius 0.1 m. Resistivity is
1e-8 ohm m. No full cell-pair matrix is retained by the integration path.

| Mesh | Current implementation | Touching analytic, direct order two | Speedup | Touching analytic, direct order three | Speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| 808 cells / 649 coefficients | 24.944 s | **4.738 s** | **5.26x** | **7.423 s** | **3.36x** |
| 1,804 cells / 1,603 coefficients | 93.135 s | **13.762 s** | **6.77x** | **27.627 s** | **3.37x** |

These are not repeated timing distributions. The policies were run current,
order two, then order three on each mesh.

## Slow rates and current patterns

All quantities below compare the experiment with the **current implementation
at rtol=1e-2 on the same mesh**, rather than an independent exact reference.
The comparison solves `R v = lambda M v` for the twenty slowest modes.

| Mesh | Direct order | Maximum relative change among first 3 rates | First 10 rates | First 20 rates |
| --- | ---: | ---: | ---: | ---: |
| 808 cells | 2 | 1.40e-5 | 3.24e-5 | **3.24e-5** |
| 808 cells | 3 | 3.15e-7 | 6.42e-7 | **6.42e-7** |
| 1,804 cells | 2 | 1.05e-5 | 2.01e-5 | **2.90e-5** |
| 1,804 cells | 3 | 2.74e-7 | 4.92e-7 | **7.23e-7** |

Current patterns are compared in the physical volume-current L2 metric, using
`R = rho * integral(J_i dot J_j) dV`. Principal angles compare subspaces and
allow rotations within nearly degenerate mode groups.

| Mesh | Direct order | Largest first-3 subspace angle | First-10 angle | First-20 angle | Minimum individual correlation among first 20 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 808 cells | 2 | 0.000992 degrees | 0.071600 degrees | **0.123420 degrees** | 0.99997620 |
| 808 cells | 3 | 0.000032 degrees | 0.000938 degrees | **0.000948 degrees** | 0.99999996 |
| 1,804 cells | 2 | 0.000619 degrees | 0.010665 degrees | **0.045373 degrees** | 0.99997340 |
| 1,804 cells | 3 | 0.000022 degrees | 0.000293 degrees | **0.001285 degrees** | 0.99999999 |

The first triplet is well separated from higher modes, with a relative gap
greater than one on both meshes. The ten-mode boundary has a smaller gap,
approximately 0.0015/0.0020 on the coarse/refined mesh. Thus changes in individual
orientations or truncated subspaces should be interpreted with their spectral
separation in mind. Twenty-one-mode boundary gaps were not measured.

## Global matrix and energy changes

Generalized eigenvalues of `M_experiment x = lambda M_current x` bound the energy
ratio over **all** coefficient vectors and test experimental positive
definiteness. All four experimental matrices remain positive definite.

| Mesh | Direct order | Relative Frobenius matrix difference | Energy ratio range over all currents | Largest relative energy change restricted to first 20 reference modes |
| --- | ---: | ---: | ---: | ---: |
| 808 cells | 2 | 2.87e-4 | 0.99885098 to 1.00110613 | 4.47e-5 |
| 808 cells | 3 | 1.21e-5 | 0.99991084 to 1.00009632 | 1.00e-6 |
| 1,804 cells | 2 | 3.49e-4 | 0.99839929 to 1.00152937 | 3.33e-5 |
| 1,804 cells | 3 | 1.54e-5 | 0.99986705 to 1.00014328 | 9.64e-7 |

The worst energy change for order two is 0.1601%; for order three it is 0.01433%.
These global differences are larger than those restricted to the smooth, slow
mode subspaces. Generalized eigenpair residuals are below 1.7e-13. Matrices are
finite and exactly symmetric.

## Contact counts and scaling

Counts are unique unordered pairs, including self interactions.

| Contact | 808 cells | 1,804 cells |
| --- | ---: | ---: |
| Self | 808 | 1,804 |
| Shared face | 1,456 | 3,406 |
| Shared edge only | 4,716 | 11,358 |
| Shared vertex only | 16,436 | 39,481 |
| All analytic contacts | **23,416 (7.16%)** | **56,049 (3.44%)** |
| Non-touching, direct quadrature | 303,420 | 1,572,061 |
| Total | 326,836 | 1,628,110 |

The current scheme uses analytic integration for 178,579/582,841 pairs, or
54.64%/35.80%. Restricting it to touching neighbors reduces that work sharply.
Analytic contact count grows 2.39x for 2.23x more cells, consistent with a local
neighbor set. Non-touching direct work remains quadratic. Shared vertices are
tested blockwise without building a dense global adjacency matrix.

## Interpretation

For the stated focus on slow, smooth patterns, fixed order two is a strong
speed candidate: 5.3–6.8x faster, first-twenty rate changes below 0.00324%, and
small current-subspace angles. Order three offers more accuracy margin while
remaining approximately 3.4x faster on both meshes.

These sphere experiments support the self/touching strategy, but do not supply
an error guarantee for arbitrary non-touching cells. Small gaps, unequal cell
sizes, distorted elements, and vertex-order dependence of low-order product
quadrature were not independently characterized here. The production distance-
controlled policy remains available as the comparison reference.

## Reproduction and verification

```sh
.venv/bin/python benchmarks/benchmark_volume_touching_policy.py
```

Optional `--orders 2 3 4` selects the fixed direct rules; `--meshes coarse` runs
only the smaller sphere. The refined archive defaults to
`benchmarks/volume_eddy_profile/refined_sphere/mesh.npz`.

The experiment also checked all five contact classes (zero through four shared
vertices) on a conforming cube mesh. On a cube whose six tetrahedra all share a
body diagonal, the experimental and current cell matrices match exactly. The
scoped selector replacement restores the original function after each run.
The JSON report stores complete settings, counts, timing observations, rates,
energy ratios, residuals, and current-pattern comparisons.
