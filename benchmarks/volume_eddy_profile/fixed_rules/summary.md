# Fixed-rule assembly cleanup

The self/touching-neighbor policy is now the production implementation:

- `quadrature_order=2`: direct double quadrature for non-touching cells.
- `touching_order=4`: fixed outer quadrature of analytic source potentials for
  self and shared-face/edge/vertex neighbors.
- No distance-based order selection, successive-order checks, refinement loop,
  `rtol`, or `max_order` arguments remain in the matrix assembly API.
- `_CellIntegrals` prepares both cell rules and source geometry once per
  streamed assembly. The point-potential stability fallback remains independent
  and prepares its source quadrature lazily.

```python
magnetic = inductance_matrix(basis)
comparison = inductance_matrix(basis, quadrature_order=3, touching_order=6)
```

The second call uses the same fixed-order code with more points. It is an
order-sensitivity check, not a separate guaranteed-accuracy assembler.

## Validation against the approved touching/adaptive experiment

Before removing the adaptive code, its direct-order-two touching policy was
assembled on both sphere meshes. Those matrices were saved locally for the
cleanup comparison. `results.json` records the resulting comparisons.

| Measurement | 808 cells | 1,804 cells |
| --- | ---: | ---: |
| Approved touching/adaptive assembly snapshot | 4.77 s | 13.64 s |
| Fixed-rule cleanup snapshot | 3.37 s | 10.08 s |
| Relative Frobenius matrix change | 3.59e-5 | 9.28e-5 |
| Maximum relative first-20 rate change | 5.12e-8 | 7.45e-8 |
| Largest first-20 current-subspace angle | 0.000226 degrees | 0.000250 degrees |
| Energy ratio range over all currents | 0.9988264 to 1.0003993 | 0.9972057 to 1.0026444 |

Both matrices remain positive definite. Removing refinement affects high
frequency current combinations more strongly than the targeted slow modes.
The first twenty sorted rates remain within 3.24e-5/2.90e-5 relatively of the
earlier distance-controlled production implementation.

A final run of the maintained refinement benchmark, after the lazy source
fallback cleanup, measured **3.40 s and 10.70 s**. Its complete settings, timings,
eigenpair residuals, contact counts, and path counters are in
`../fixed_rules_refinement/results.json`. These are single-run observations.

## Checks and maintained tools

- All **114 volume tests pass**, including independent point/self/cube integral
  references, conforming contact dispatch, reciprocity, energy, scaling, and
  generalized modes. Strict reference tests request fixed higher orders.
- CPU, line, memory, and baseline-comparison benchmark smoke runs passed.
- The fixed-order mode-comparison runner passed on the 808-cell sphere against
  direct order three/touching order six, with positive energy ratios and
  eigenpair residuals below 4.3e-14.
- Python compilation, formatting, and diff checks passed.

Current runners:

```sh
.venv/bin/python benchmarks/benchmark_volume_sphere_refinement.py
.venv/bin/python benchmarks/compare_volume_eddy_orders.py
.venv/bin/python benchmarks/profile_volume_eddy_assembly.py
```

Examples now expose `--quadrature-order` and `--touching-order` where applicable;
`--check-convergence` repeats assembly at higher fixed orders. Historical
adaptive reports are retained, while their obsolete comparison runners have
been removed. An accurate matrix assembler is deferred to future work.
