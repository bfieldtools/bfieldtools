# Volume implementation and validation

This audit maps the current implementation to its tests and examples. Test
coverage below describes assertions in the repository, rather than a claim
that a particular numerical result is a physical error bound. The remaining
work is summarized in [PLAN.md](PLAN.md); the equations and magnetization
interpretation are in [THEORY.md](THEORY.md).

The five volume test modules passed together during this documentation audit:
99 tests passed, including the optional TetGen sphere test.

## Evidence by implementation area

| Implementation | Evidence | What it establishes |
| --- | --- | --- |
| [`mesh.py`](mesh.py) | [`test_volume_mesh.py`](../../tests/test_volume_mesh.py) | Orientation, volumes, quality, local incidence, invalid input/backend handling; optional TetGen sphere boundary/volume preservation and refinement |
| [`topology.py`](topology.py) | [`test_volume_topology.py`](../../tests/test_volume_topology.py) | Incidence identities, physical face-flux signs, zero boundary flux, independent and complete cotree currents against the conserved face-flux space; disconnected solid balls; rejection of handles/cavities/nonmanifold boundaries |
| [`basis.py`](basis.py), [`impedance.py`](impedance.py) resistance | [`test_volume_resistance.py`](../../tests/test_volume_resistance.py) | Analytic cube-loop resistance; Joule loss and bilinear forms against independent barycentric reconstruction; current continuity; material and length scaling; complex/multiple coefficient fields; basis/layout contracts |
| [`integrals.py`](integrals.py), [`impedance.py`](impedance.py) inductance | [`test_volume_inductance.py`](../../tests/test_volume_inductance.py) | Independent point/self/pair integral references; touching/near/far cells; rigid motions and scaling; blocked versus dense contraction; energy; basis-invariant generalized modes and metric checks; fixed-order contact dispatch and independent high-order references |
| [`dynamics.py`](dynamics.py) | [`test_volume_dynamics.py`](../../tests/test_volume_dynamics.py) | Independent wire-kernel/coupling integration, orientation and intersection rejection, analytic ramp/free-decay solutions, modal energy balance, and bundled-mesh contract |

The topology tests compare the retained curls with the full nullspace of the
interior-face divergence matrix. They do not construct an explicit dual-cycle
basis. Boundary-edge removal enforces zero tangential trace of the potential;
existing numerical assertions check the resulting current fluxes. Direct tests
of potential traces, total dipole moment, and magnetization-based coupling
would accompany future potential/dipole operators.

The general `CellCurrentBasis` validates operator storage and shape. It accepts
operators supplied by callers, so its container checks alone establish neither
current conservation nor independence.

## Reproducible examples

Run from the repository root, with optional meshing dependencies installed:

```sh
uv run --locked --extra meshing python examples/volume/sphere_meshing.py
uv run --locked --extra meshing python examples/volume/topology.py
uv run --locked --extra meshing python examples/volume/resistance.py
uv run --locked --extra meshing python examples/volume/eddy_modes.py
uv run --locked --extra meshing python examples/volume/eddy_currents.py
```

`sphere_meshing.py` demonstrates surface-to-volume generation and quality
inspection. The remaining examples reuse the bundled unit-radius sphere:
226 vertices, 808 tetrahedra, and 649 independent current coefficients, loaded
with `load_example_tetra_mesh()`. Scale vertex coordinates to physical metres.
Refining this volume mesh alone does not remove the faceted-surface geometry
error.

`topology.py` displays boundary, tree, and cotree edges plus a selected edge's
circulating current. `resistance.py` displays ordinary resistance eigenpairs;
their Euclidean coefficient normalization makes them gauge-dependent
diagnostics. Physical decay modes require the generalized inductance metric,
as demonstrated by `eddy_modes.py`.

`eddy_currents.py` uses a sphere of radius 0.1 m with resistivity 1e-8 Ω m,
driven by an external loop of radius 0.075 m at z=0.14 m. It applies a smooth
1000 A pulse lasting 80 ms and displays the waveform, decay times, modal
amplitudes, and current snapshots at 20, 60, and 180 ms. All modes are retained;
this example does not measure modal-truncation error or solve source-circuit
back reaction.

On each run, the transient example checks generalized eigenpairs,
M-orthogonality, positive rates, coupling quadrature refinement, physical
face-flux continuity/boundary leakage, free decay and decreasing energy, and
agreement with direct integration of the same coefficient ODE. Its
`--check-convergence` option also increases the two fixed quadrature orders and
doubles the loop segmentation. `eddy_modes.py` has a quadrature refinement option
with the same name. Neither performs spatial mesh refinement.

The topology, resistance, mode, and transient examples accept
`--screenshot /tmp/volume-example.png` for offscreen output. Resistance and
mode examples accept `--modes` to select three indices. Use each script's
`--help` for its controls.

## Acceptance for physical validation

1. Compare a sphere or slab against an independent analytical decay/driven
   reference and a thin conductor against the existing surface formulation.
2. Report low-mode rates and driven observables under both boundary-geometry
   and volume-mesh refinement; resolve skin depth in the tested regime.
   Compare invariant subspaces where symmetry makes eigenvalues degenerate.
3. Refine quadrature, loop segmentation, time integration, and retained mode
   count independently so each source of error is visible. The current
   all-mode/direct comparison verifies evolution of the discrete system.
4. Validate new field operators at exterior, near-boundary, and interior points
   against independent current integration, including difficult thin cells.
5. Preserve existing conservation, boundary, reciprocity, scaling, energy,
   basis-invariance, and eigenpair checks throughout these studies.

For executable regression checks:

```sh
uv run --locked --extra dev --extra meshing pytest tests/test_volume_mesh.py tests/test_volume_topology.py tests/test_volume_resistance.py tests/test_volume_inductance.py tests/test_volume_dynamics.py
```
