# Reducing integral array copies

Compared with commit `c163220`, with unchanged chunk sizes and quadrature settings:

- Make tetrahedron face geometry contiguous once, so the displacement reshape
  shares memory instead of copying each point/source block.
- In `gamma0` and `x_distance`, contract displacement vectors first, then roll
  scalar results. Only small geometry arrays are rolled in coordinate form.
- Use a vertex slice in `d_distance` instead of a copying `np.take`.

Five warmed, alternating-order benchmark repetitions:

| Full inductance assembly | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| 48 cells / 25 DOFs | 345.87 ms | 328.05 ms | 1.05x |
| 162 cells / 109 DOFs | 3.939 s | 3.737 s | 1.05x |

All nine kernel/interaction cases and both complete inductance cases agreed
exactly with the baseline. The comparison loads both volume and shared triangle
code from the baseline revision.

Separate default CPU profiles: reshape self time drops from 0.289 s to 0.0033 s;
NumPy roll cumulative time drops from 0.912 s to 0.602 s. Profiling overhead is
included in these function times, so use the benchmark above for overall speed.
Unprofiled process-CPU medians in the separate profiling runs are 3.902 s before
and 3.656 s after, with normal timing variation between runs.

Peak traced allocations remain about 5.7 MiB for the default 162-cell assembly.
The change reduces transient copies and data movement but does not lower the
maximum live allocation reached by other stages of the calculation. The larger
64/1024/64 block/point/source configuration still peaks at about 86.6 MiB.

Validation: 113 integral, magnetic-field, and volume tests passed. The edge
integral now has an independent numerical quadrature check on read-only,
noncontiguous inputs, with symmetrization enabled and disabled. A full-suite
attempt terminated with a VTK segmentation fault while rendering in
`test_line_conductor_functionality`; the numerical subset completed normally.

Artifacts: `volume_copy_reduction_results.json`,
`volume_copy_reduction_memory_results.json`, and `volume_cpu_copy_reduction/`.
