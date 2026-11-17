# Surface integral benchmarks

Functional versus cached TriangleIntegrals on flat triangle meshes and a curved triangle surface.
Timings are medians of five warmed calls with alternating execution order.
The cached backend is substituted only inside this benchmark; production surface code is unchanged.

Full self-inductance includes quadrature, near/far classification where enabled, analytic diagonal correction, and projection to vertex stream functions. Mesh preparation is excluded.
Assembly uses the six-point Dunavant rule (rule index 3), point chunks of approximately 128, margin 2, and deterministic non-clustered chunks.

| Workload | Faces | Functional ms | Cached ms | Speedup |
| --- | ---: | ---: | ---: | ---: |
| disc: uniform kernel (planar) | 143 | 2.20 | 2.14 | 1.03x |
| disc: uniform kernel (off-plane) | 143 | 3.58 | 3.35 | 1.07x |
| disc: linear dipole kernel | 143 | 4.16 | 3.85 | 1.08x |
| disc: self-inductance (exact) | 143 | 10.34 | 10.20 | 1.01x |
| disc: self-inductance (approx far) | 143 | 17.46 | 17.25 | 1.01x |
| plate: uniform kernel (planar) | 1258 | 19.53 | 19.52 | 1.00x |
| plate: uniform kernel (off-plane) | 1258 | 30.69 | 28.19 | 1.09x |
| plate: linear dipole kernel | 1258 | 34.96 | 32.93 | 1.06x |
| plate: self-inductance (exact) | 1258 | 672.10 | 667.27 | 1.01x |
| plate: self-inductance (approx far) | 1258 | 390.90 | 388.68 | 1.01x |
| sphere: uniform kernel (curved) | 320 | 7.61 | 6.99 | 1.09x |
| sphere: linear dipole kernel | 320 | 8.61 | 8.20 | 1.05x |
| sphere: self-inductance (exact) | 320 | 69.06 | 64.20 | 1.08x |
| sphere: self-inductance (approx far) | 320 | 57.82 | 55.23 | 1.05x |

The planar uniform kernel omits the solid-angle term, so it cannot benefit from sharing vertex distances between gamma0 and omega. Its timings are essentially unchanged. Off-plane and curved uniform evaluations improve by about 7–9%; linear dipole kernels improve by about 5–8%.

All compared uniform kernels and self-inductance matrices are identical. The maximum relative Frobenius difference in the dipole kernels is about 1.1e-15. Both backends use the same approximation and analytic diagonal settings.

Reproduce from the repository root:

```sh
MPLCONFIGDIR=/private/tmp/bfieldtools-mpl .venv/bin/python benchmarks/benchmark_surface_integrals.py --output benchmarks/surface_integrals_results.json
```

Raw timing samples and numerical differences are in surface_integrals_results.json.
