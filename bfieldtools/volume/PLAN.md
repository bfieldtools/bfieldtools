# Volumetric eddy-current modeling: status and plan

The first loop-driven sphere demonstration is implemented: mesh → independent
volume currents → resistance and inductance → decay modes → driven transients.
Spatial convergence and independent physical validation remain the next milestone.
This status is based on the source, tests, and examples listed below; an example's
quadrature checks do not establish spatial accuracy.

## Documentation map

- [THEORY.md](THEORY.md): equivalent magnetization and the existing surface
  dipole-layer formulation, edge loops, compatible spaces, topology, and dynamics.
- [EDGE_CURRENT_MATRICES.md](EDGE_CURRENT_MATRICES.md): local normalization,
  opposing-edge curls, resistance/inductance assembly, and tree–cotree gauge.
- [INDUCTANCE_INTEGRALS.md](INDUCTANCE_INTEGRALS.md): tetrahedron potentials,
  implemented numerical integration, accuracy limits, and integration research.
- [ARCHITECTURE.md](ARCHITECTURE.md): containers, dependency boundaries, and APIs.
- [VALIDATION.md](VALIDATION.md): implementation evidence, runnable examples,
  existing checks, and acceptance criteria for the next milestone.

## Implemented scope

Isolated stationary conductors, no terminals, scalar positive resistivity uniform
or constant per tetrahedron, permeability μ₀, and magnetoquasistatic dynamics.
The current approximation is constant per tetrahedron with matching normal
components and zero boundary leakage. The equivalent magnetization **T** is a
lowest-order Nédélec field with zero tangential boundary trace; its gauge is
fixed by the tree–cotree construction.

| Area | Implemented | Current limit |
| --- | --- | --- |
| Geometry (`mesh.py`) | Validated array container, oriented cells/boundary, volumes, quality, PyVista conversion, bundled sphere | Input arrays must already be geometrically conforming; no overlap/self-intersection/hanging-node detection |
| Meshing (`mesh.py`) | Optional TetGen adapter for one closed outward genus-zero surface | No cavities, multiple solids, labels, material interfaces, or local size fields; maximum volume is a target |
| Topology (`topology.py`) | Incidence operators and independent boundary-compatible tree–cotree gauge | Gauge supports disconnected solid-ball components; rejects handles, cavities, and nonmanifold boundaries |
| Current representation (`basis.py`) | General cell-current operators, edge-curl builder, real/complex and multiple-field reconstruction | No public potential/dipole-density reconstruction or explicit dual-cycle/global-pattern builder |
| Resistance (`impedance.py`) | Exact sparse weighted Gram assembly | Scalar resistivity only |
| Inductance (`integrals.py`, `impedance.py`) | Near-source analytic potentials, numerical outer integration, reciprocal cell blocks, dense basis matrix, optional cached cell interactions | Empirical quadrature estimator; no general exact pair evaluator, spatial adaptivity, compression, or matrix-free solve |
| Modes (examples using SciPy) | Dense generalized symmetric eigensolves, residual and metric checks | No library eigensolver helper or large-system solver |
| Excitation/evolution (`dynamics.py`) | External polygonal-loop coupling and prescribed-drive modal integration | No general excitation API or coupled source-circuit dynamics; wires must remain outside the conductor |

The distinction between the meshing adapter and topology support matters:
multiple supplied solid-ball meshes already work with the gauge, even though the
adapter cannot generate multiple solids in one call. Conservation and boundary
conditions are guaranteed by the edge construction, not by arbitrary operators
passed to `CellCurrentBasis`.

## Next milestone: physical accuracy

1. **Spatial convergence and independent references.** Refine volume and surface
   geometry, compare sphere/slab decay or driven response with analytical
   solutions, and compare a thin conductor against bfieldtools' surface model.
   Track low eigenspaces for degenerate modes. Separate spatial, geometry,
   quadrature, loop-segmentation, time-integration, and modal-truncation errors.
   Include skin-depth resolution and a mode-count convergence study.
2. **Volume-current field evaluation.** Add vector-potential and magnetic-field
   operators, tested outside, near, and inside cells against independent
   integration. The scalar `tetrahedron_potential` kernel is available, but a
   general magnetic-field operator is not. Add potential reconstruction and a
   total-dipole-moment operator if using the magnetization interpretation for
   diagnostics or alternative coupling.
3. **Difficult integration cases.** Validate thin cells and extreme scales;
   introduce spatial subdivision when order refinement fails or is inefficient.
   Expand independent high-accuracy references and control far-source accuracy
   separately. See the integration document for the detailed backlog.

Acceptance criteria and current evidence are in [VALIDATION.md](VALIDATION.md).

## Subsequent extensions

- **General topology:** cavity-compatible gauge constraints, handle-spanning
  global currents, and an explicit face-flux/dual-cycle reference builder.
  Validate completeness and independence, beyond the current solid-ball tests.
- **Geometry workflows:** robust intersection/overlap/conformity validation,
  cavities, multiple-region meshing, labels, and local refinement. Add block,
  cavity, and torus examples as support becomes available. TetGen is the current
  backend; additional adapters are optional follow-up work.
- **Materials and physics:** tensor resistivity, magnetic materials, terminals
  and imposed transport currents, source-circuit back reaction, and higher-order
  elements. These require extensions to the present model.
- **Scale and integration:** profile assembly/eigensolves after validation;
  investigate compressed or matrix-free inductance and partial eigensolvers.
  Add a higher-level conductor API and public documentation once interfaces
  stabilize. `magnetics.py` and `modes.py` remain possible future modules;
  `meshing.py` does not exist separately from `mesh.py`.
