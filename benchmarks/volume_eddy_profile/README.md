# Volume assembly experiments

The current implementation uses fixed direct order two and analytic
self/shared-vertex order four. It has no tolerance-driven refinement.

Current runners:

- `benchmarks/benchmark_volume_sphere_refinement.py`: coarse/refined mesh timings
  and path counts with the fixed rules.
- `benchmarks/compare_volume_eddy_orders.py`: compare fixed orders two/four with
  three/six, including slow rates, current subspaces, and magnetic energy.
- `benchmarks/profile_volume_eddy_assembly.py`: fixed-rule sphere CPU profile.

Older JSON/text/profile outputs and their summaries document the preceding
adaptive and distance-bound experiments. Their settings and reproduction
commands describe that earlier source version. The obsolete tolerance-comparison
and patched-policy runners were removed when the fixed policy was adopted;
the recorded results are retained as historical evidence.

`refined_sphere/mesh.npz` is the archived 1,804-cell refinement and is reused by
the current runners. `fixed_rules/results.json` records the final simplification
check against pre-change touching/adaptive matrices saved during the cleanup.
