# Triangle integration API: decision deferred

Status: retain both implementations for now; revisit API consolidation later.

The functional interface in `integrals.py` remains the surface code's API.
`integrals_cached.py` contains the `TriangleIntegrals` dataclass, derived from
those formulas. Volume integration uses the dataclass. Its cached properties
share geometry, vertex distances, and face-plane distances within an evaluation
batch. Optional supplied edges and areas seed their respective caches. Instances
borrow their input arrays and are discarded after each volume evaluation chunk.

This split is intentional at this checkpoint. We have not decided whether the
cached class should become the common triangle-integration API, remain a volume
implementation detail, or back compatible functional wrappers. Surface callers
have not been migrated as part of this work.

## Evidence at this checkpoint

- The cached volume implementation reduced a 162-cell / 109-DOF full inductance
  assembly from 3.75 s to 3.46 s, with identical benchmark outputs.
- Peak traced allocations for the default assembly were approximately 5.9 MiB,
  versus the earlier functional measurement of approximately 5.7 MiB.
- Planar surface evaluations gain little because they omit the solid-angle
  term. Off-plane and curved-surface kernels benefit from shared distances;
  measured curved-surface assembly speedup was about 8% for the tested mesh.
- Numerical checks cover the functional and cached exact, planar, normalized,
  approximate, dipole, and vertex formulas, as well as supplied read-only edges,
  noncontiguous inputs, regularization options, and cache reuse.

See `../benchmarks/volume_cached_integrals_summary.md` and
`../benchmarks/surface_integrals_summary.md` for configurations and measurements.
These results describe the measured workloads, not a general speed guarantee.

## Questions for the later API decision

1. Should the public interface expose cached properties, functions, or both?
2. Should source-only geometry have a separate lifetime from evaluation-point
   data, allowing reuse across batches without retaining displacement arrays?
3. How should supplied edges/areas, planar evaluation, regularization, and
   normalized versus unnormalized projections be represented consistently?
4. What ownership and mutation contract should apply to borrowed input arrays
   and cached results, and when should surface callers migrate?

Until that decision, preserve the functional API and keep numerical comparisons
between the two implementations. Further formula fusion and algebraic changes
can be evaluated separately from API consolidation.
