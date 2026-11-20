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
| [`mesh.py`](mesh.py) | [`test_volume_mesh.py`](../../tests/test_volume_mesh.py) | Orientation, volumes, quality, local incidence, invalid input/backend handling; optional TetGen sphere/torus boundary and volume preservation and refinement; torus tunnel and global-mode count; rejection of disconnected/cavity boundaries |
| [`topology.py`](topology.py) | [`test_volume_topology.py`](../../tests/test_volume_topology.py) | Incidence identities, physical face-flux signs, zero boundary flux, independent and complete cotree currents against the conserved face-flux space; disconnected solid balls; rejection of handles/cavities/nonmanifold boundaries |
| [`basis.py`](basis.py), [`impedance.py`](impedance.py) resistance | [`test_volume_resistance.py`](../../tests/test_volume_resistance.py) | Analytic cube-loop resistance; Joule loss and bilinear forms against independent barycentric reconstruction; current continuity; material and length scaling; complex/multiple coefficient fields; basis/layout contracts |
| [`integrals.py`](integrals.py), [`impedance.py`](impedance.py) inductance | [`test_volume_inductance.py`](../../tests/test_volume_inductance.py) | Independent point/self/pair integral references; touching/near/far cells; rigid motions and scaling; blocked versus dense contraction; energy; basis-invariant generalized modes and metric checks; fixed-order contact dispatch and independent high-order references |
| [`dynamics.py`](dynamics.py) | [`test_volume_dynamics.py`](../../tests/test_volume_dynamics.py) | Independent wire-kernel/coupling integration, orientation and intersection rejection, analytic ramp/free-decay solutions, modal energy balance, and bundled-mesh contract |

The topology tests compare the retained curls with the full nullspace of the
interior-face divergence matrix. They use an independent dense divergence-nullspace
reference. Boundary-edge removal enforces zero tangential trace of the potential;
existing numerical assertions check the resulting current fluxes. Direct tests
of potential traces, total dipole moment, and magnetization-based coupling
would accompany future potential/dipole operators.

The general `CellCurrentBasis` validates operator storage and shape. It accepts
operators supplied by callers, so its container checks alone establish neither
current conservation nor independence.

## Complete current spaces on handles and cavities

`test_volume_currents.py` checks the stacked-nullspace augmented
construction on solid balls, a solid ring, an enclosed cavity, a two-handle
solid, a solid with both a handle and a cavity, disconnected rings with an unused
vertex, a thin two-handle slab with a wide stack, and a single tetrahedron. It
checks independence and completeness against the full conserved-flux space,
and integrates the
reconstructed Cartesian current through every cell face. The ring's global
mode has nonzero cross-sectional cut flux; all local curl columns have zero
net cut flux. A cavity adds no handle mode.

Resistance, inductance, generalized decay rates, and external-loop coupling
agree under change of basis against an independently constructed full
divergence-nullspace ring basis. Tests cover complex coefficients, energy, permutations/orientations, geometric scaling,
nonmanifold edge/vertex rejection, invalid supplied fluxes, and the dense size
guard. These establish algebraic completeness and assembly compatibility, not
spatial convergence of a torus response. The original tree gauge still rejects
handles and cavities by design; callers use the augmented constructor instead.

Run the backend-independent example with
`uv run --locked python examples/volume/global_currents.py`; append `--show`
to display the ring current. Run the new checks with
`uv run --locked --extra dev pytest tests/test_volume_currents.py`.

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

## Meshed torus transient

Run the optional-TetGen example with

```sh
uv run --locked --extra meshing python examples/volume/torus_eddy_currents.py
```

Append `--no-show` for numerical checks, `--screenshot /tmp/torus.png` to export
the figure, or `--check-convergence` to report changes with tighter inductance
quadrature and refined source-loop integration. `--major-sections`,
`--minor-sections`, and `--max-volume` control the geometry/volume mesh.

The example meshes a torus with major radius 0.1 m and tube radius 0.025 m.
Resistivity is 1e-8 ohm metres; a coaxial polygonal loop of radius 0.15 m at
z=0.06 m carries a prescribed 1000 A, 40 ms smooth pulse. It uses all modes
from the augmented basis, verifies exactly one global current degree of
freedom, and reconstructs the driven currents and subsequent free decay.

On each run it checks generalized eigenpairs, M-orthogonality, modal versus
direct coefficient-space integration, decreasing magnetic energy after the
pulse, and physical face-flux reconstruction. Integrating cell currents through
two oblique tube cuts produces the same circulating current; local curl
columns have zero net cut flux, while the drive excites the global current.
The PyVista window embeds the source current, circulating current, and magnetic
energy traces above three linked three-dimensional snapshots. Each snapshot
removes the upper front quarter to expose surface and interior cell magnitudes,
using one shared absolute color scale. Cyan arrows sample the original current
at uniformly spaced physical locations around the tube. Each vector comes from
its containing tetrahedron, with the same sampling positions in all snapshots.
Arrows have fixed length to emphasize direction; the cell colors carry current
magnitude.
`--show-mesh` reveals cell edges; the default keeps them hidden for readability.
The coarse default mesh and
integration sensitivity checks do not establish spatial/skin-depth convergence
or an independently validated torus response.
