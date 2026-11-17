# CPU profile of streamed inductance assembly

Workload: unit cube split into 162 tetrahedra, 109 current DOFs. Mesh and basis
construction are excluded. Integration settings: order 4, rtol 1e-4, max_order 12.
The default block/source/point limits are 16/16/256; the alternative uses
64/64/1024. Results agree within 2.2e-15 relative Frobenius error.

Unprofiled medians of three warm calls:

| Configuration | Process CPU | Wall time |
| --- | ---: | ---: |
| Default | 3.902 s | 3.904 s |
| Larger batches | 3.815 s | 3.818 s |

Separate cProfile runs use `time.process_time`; profiler overhead increases the
default run to 4.422 s and the larger-batch run to 3.865 s. Use the unprofiled
medians for performance comparisons. Percentages below refer to the default
profile, and the three rows do not overlap.

| Stage | Profiled CPU | Share |
| --- | ---: | ---: |
| Triangle potential | 3.471 s | 78.5% |
| Other geometry integration | 0.852 s | 19.3% |
| Projection into basis and other assembly overhead | 0.099 s | 2.2% |

Within the triangle kernel, cumulative `gamma0` time is 1.732 s (39.2% of full
assembly), `omega` is 1.014 s (22.9%), `x_distance` is 0.483 s (10.9%), and
`d_distance` is 0.146 s (3.3%). NumPy `roll` takes 0.912 s cumulatively across
18,117 calls (20.6% of assembly); it is already included in those functions.
Native `c_einsum` takes 1.509 s across 24,882 calls, also included above.

The independent order instrumentation found only mandatory order-4 and order-6
integration: 121 calls at each order for the default, nine at each order for
larger batches. Every directed cell pair is evaluated at both orders. Order 4
evaluates 1,679,616 point/source pairs and takes about 0.897 CPU seconds; order 6
evaluates 5,668,704 pairs and takes about 3.002 CPU seconds. No order-8 or higher
adaptive refinement occurs for this mesh. Adaptive-row batching is therefore
not a useful target for this benchmark. The workload does not exercise the
far-source quadrature branch.

Larger batches reduce triangle calls from 2,013 to 141, but runtime falls only
2.2%. The profile points to array computations and copies inside the kernel,
rather than call overhead, as the next optimization targets. It does not by
itself distinguish arithmetic throughput from memory-bandwidth limits.

An additional array-layout check found that `vertices[:, TETRA_FACES]` produces
noncontiguous faces. With 16 sources, face strides are `(24, 1152, 384, 8)`;
subtraction produces noncontiguous displacements. Flattening source and face
axes consequently copies the displacements (`np.shares_memory` is false).
Making the faces contiguous before subtraction makes that reshape a view.
Reshapes called from `tetrahedron_potential` account for 0.288 CPU seconds
(6.5% of assembly). This is a concrete candidate to benchmark; no kernel change
has been made here. Other candidates are sharing vertex distances between
`gamma0` and `omega`, sharing edge geometry, and reducing full-array `roll`s.

Reproduce from the repository root:

```sh
MPLCONFIGDIR=/private/tmp/bfieldtools-mpl .venv/bin/python benchmarks/profile_volume_cpu.py
```

`default.txt` and `larger_batches.txt` contain cumulative and self-time tables;
`.prof` files can be opened with pstats or a compatible viewer. `results.json`
contains timings, full function statistics, and quadrature counters.
