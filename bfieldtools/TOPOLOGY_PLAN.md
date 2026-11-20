# Plan: complete volume and surface current spaces

Status: volume reference methods are implemented; surface support and the sparse
local/global decomposition remain planned. This plan
extends the volume [roadmap](volume/PLAN.md) and the existing surface
stream-function construction. Both cases must support handles, even though
closed surfaces such as a torus are uncommon in present surface workflows.

Implemented volume APIs: `augmented_flux_basis` supplies a size-limited dense
stacked-nullspace decomposition into independent interior edge curls and global
modes; `build_face_current_basis` reconstructs the fluxes for impedance and
dynamics assembly. Cavity gauge redundancy is removed by rank selection. The TetGen adapter accepts handles and the torus
eddy-current example includes excitation and decay. A general cut-normalization
API, sparse augmented construction, and cavity meshing remain planned. See
[ARCHITECTURE.md](volume/ARCHITECTURE.md#handles-and-cavities) for usage and
[VALIDATION.md](volume/VALIDATION.md#complete-current-spaces-on-handles-and-cavities)
for current evidence. The delivery sequence below remains the design target;
it does not imply all stages are complete.

## Scope and conventions

Support stationary isolated conductors with no terminals: currents conserve
charge locally and cannot cross a conductor boundary. Accept orientable,
conforming manifold meshes, including disconnected components. Reject invalid
topology explicitly. Keep geometric intersection validation a separate concern.

Distinguish three cases:

- A hole in an open surface introduces an additional boundary loop. The current
  surface `inner` basis already gives each such loop one constant stream-function
  value, with one reference boundary fixed per component.
- A handle introduces global circulation that local potential derivatives can
  miss: one additional volume-current mode for a solid torus, and two additional
  surface-current modes for its closed boundary surface.
- An enclosed volume cavity requires extra potential gauge constraints, but
  does not by itself introduce a handle current. A spherical shell has no missing
  global current mode of this type.

Use oriented integer incidence matrices, separate from metric-dependent FEM
operators:

| Domain | G | C | D |
| --- | --- | --- | --- |
| Volume | vertex-to-edge differences | edge-to-face circulation | face-to-cell balance |
| Surface | vertex-to-edge differences | edge-to-triangle circulation | absent |

Require CG = 0 and, in the volume, DC = 0. These are not automatically the
Cartesian component matrices exposed by `mesh_calculus.gradient_matrix` or
`curl_matrix`; introduce explicit incidence builders with documented signs.

## Shared numerical reference

For edge cochains, ker(C) / im(G) identifies curl-free fields that are not
single-valued gradients. Euclidean representatives satisfy CH = 0 and G^T H = 0.
This detects absolute first cohomology; its dimension is one in a solid torus
and two on a closed torus surface. It does not directly construct volume face
currents, and boundary conditions change the relevant complex.

Implement small-mesh references using SVD or rank-revealing QR. Report matrix
shapes, rank thresholds, residuals, and inferred dimensions. Derive expected
dimensions independently from connectivity, Euler characteristics, and graph
cycle counts. Avoid forming dense nullspaces in the production path.

For a conserved space ker(A) and an existing local-current map L satisfying
AL = 0, extract a complementary basis H from

    AH = 0,  L^T H = 0.

Then [L_independent, H] must be independent and span ker(A). A positive-definite
metric W can instead select representatives through L^T W H = 0. Topology fixes
the dimension, while the metric changes the representatives. Keep topological
generators distinct from metric-dependent harmonic representatives.

## Volume construction

### Complete face-flux reference

1. Fix all boundary face fluxes to zero, including cavity walls. Let D_i contain
   the interior-face columns of D.
2. For independent small-mesh validation, construct Z spanning ker(D_i) using
   a dense divergence-only nullspace solve.
3. Verify the dimension N_interior_faces - N_cells + N_components.
4. Reconstruct the fluxes with lowest-order Raviart-Thomas functions and expose
   the resulting Cartesian cell-current operators as `CellCurrentBasis`.
   Conservation cancels the linear RT0 term, leaving constant current per cell.

The implemented augmented builder below constructs the local/global space
directly from the sparse divergence/transpose-curl stack.

### Local edge curls plus global currents

Let C_i be C restricted to interior faces and interior edges. Boundary edge
potential coefficients remain zero. It satisfies D_i C_i = 0.

1. Select independent columns of C_i, initially with rank-revealing algebra.
   Extend the tree-cotree construction to account for cavity-induced kernel
   modes; removing only gradients of zero-boundary nodal potentials is not
   sufficient. Potentials constant at different values on different boundary
   components can produce additional curl-free relative edge fields.
2. Extract H_v representing ker(D_i) / im(C_i). For a Euclidean reference solve
   the right nullspace of the sparse stack [D_i; C_i^T]. The current reference
   uses dense SVD with a size guard; sparse nullspace extraction remains planned.
   This enforces D_i H_v = 0 and C_i^T H_v = 0 with the identity flux metric.
   RT0 mass weighting can later replace C_i^T by C_i^T W.
   This is relative second cohomology, whose
   dimension is b1 of an orientable three-dimensional conductor.
3. Form q_i = C_independent a + H_v h and reconstruct cell currents using the
   same face-flux reconstruction as the complete reference.
4. Normalize global coefficients by signed flux through independent
   cross-sectional cuts. Verify the cut-flux matrix is nonsingular before
   normalization. Document orientation and units (amperes).

A solid torus must contribute one independent handle current. A spherical
cavity must change the gauge treatment without adding a handle current. Never
enable general topology merely by removing the current rejection checks.

## Surface construction

### Preserve existing boundary-hole behavior

Keep `inner2vert` and the current `inner` basis semantics. Each interior vertex
has a stream-function coefficient, each non-reference boundary loop has one
constant coefficient, and the reference boundary is zero. For closed components,
remove one constant stream-function gauge when an independent basis is needed.
The unrestricted `vertex` basis retains its current boundary semantics.

### Conserved edge fluxes and the missing handle modes

Use oriented edge-integrated conormal current fluxes q. Choose the conormal
orientation consistently so that the rotated gradient

    K = grad_s(psi) cross n

maps to q = G psi (or a documented global sign). Triangle current balance is
then Cq = 0. Derive and test this orientation relation against the existing
rotated-gradient implementation rather than assuming its sign.

1. Set flux through every open boundary edge to zero. Let C_i retain the
   interior-edge columns of C. Construct the complete reference ker(C_i)
   using an independent dense nullspace solve on small meshes.
2. Let S map existing boundary-compatible stream coefficients to vertex values.
   Define L_s = G[interior_edges, :] S. Remove constant gauges and any redundant
   columns to obtain an independent stream-current map.
3. Extract H_s representing ker(C_i) / im(L_s), initially by
   C_i H_s = 0 and L_s^T H_s = 0. Compare with the full reference nullspace.
4. Reconstruct each flux basis function as a tangent, piecewise-constant current
   per triangle, using surface RT0 reconstruction. Check conormal continuity
   across shared edges and zero flux across open boundaries.
5. Normalize the added coefficients using independent cut fluxes or equivalent
   periods of the rotated current one-form. Document which quantity is used.

For a connected orientable surface of genus g, this adds 2g modes to the
boundary-compatible stream-function space. With b > 0 boundary loops, that
stream space already includes b - 1 boundary-hole circulation coefficients.
Do not append all relative cohomology modes to it: that would count those hole
currents twice. A closed torus adds two modes, an annulus adds none, and a
punctured torus still adds two.

### Surface API and operator integration

Introduce an opt-in current basis that combines the existing stream-current
operators with global columns. Keep legacy `StreamFunction` coefficients and
ordinary genus-zero workflows compatible. Store global coefficients separately
or expose explicit coefficient slices; they are not ordinary vertex values.

Audit resistance, inductance, magnetic/vector-potential coupling, excitation,
optimization, and current reconstruction for assumptions that every current is
a rotated gradient. Add operator assembly accepting general tangent triangle
currents, then confirm agreement with existing operators on the stream subspace.
Retain local/global cross terms in all assembled matrices.

The surface `suh` basis should have an explicit policy: retain the requested
number of stream modes and append the complete handle subspace. Verify boundary
hole modes are retained or represented under truncation rather than silently
claiming completeness for a truncated basis.

Scalar dipole-layer potentials and stream-function contour extraction require
special treatment. Global currents may require multivalued potentials or cuts;
they cannot be silently passed through existing scalar stream-function APIs.
Provide direct current-based field evaluation first. Specify cut/jump-aware
potential and contour reconstruction separately before exposing them as supported.
The same restriction applies to volume magnetization/potential diagnostics.

## Sparse implementation after the reference

Build reusable connectivity and boundary-component utilities for both
dimensions. Keep the stacked operators sparse and compute only the small
nullspace, with residual and dimension checks and explicit rank tolerances.
Use sparse cavity-compatible selection of independent local curls. Sparse
constrained metric solves can improve representatives when needed.

Keep connectivity and incidence operators geometry-independent. Reconstruct
Cartesian currents and apply mass/resistance metrics in a separate layer.
Retain generator/cut metadata for interpretation without making impedance or
dynamics assembly depend on a topology object. Reuse `CellCurrentBasis` for
volume integration; design the corresponding surface adapter around existing
surface callers before committing a public API.

## Delivery sequence and acceptance

1. **Incidence and fixtures:** surface incidence builders, boundary/component
   metadata, and small deterministic fixtures: disk, annulus, sphere surface,
   closed torus surface, punctured torus, solid ball, solid ring, and spherical
   cavity. Include disconnected combinations. Fixtures should not require an
   optional meshing backend.
2. **Complete references:** independent divergence nullspaces, stacked-nullspace
   extraction, and conservative surface/volume reconstruction. Record ranks and
   expected global counts; compare represented spaces rather than exact columns.
3. **Volume support:** independent local curls, cavity gauge treatment, global
   handle columns, and existing resistance/inductance/coupling integration.
   Verify agreement with the full divergence nullspace under change of basis.
4. **Surface support:** opt-in augmented current basis, general-current operator
   assembly, and explicit behavior for `suh`, scalar potentials, and contours.
   Demonstrate the two global torus currents alongside the ordinary annulus.
5. **Sparse construction and geometry:** replace dense construction in the
   production path, measure memory/conditioning, and extend the TetGen adapter
   separately for handles and cavity seeds. Do not claim meshing support solely
   because supplied tetrahedral meshes now work.

At each stage verify incidence identities, conservation, boundary leakage,
normal/conormal continuity, independence, completeness, and cut-flux rank.
Check invariance under index permutations, orientation changes, and rigid
motions. Preserve genus-zero results and boundary-hole behavior. Compare energy,
coupling, and generalized decay modes between equivalent complete bases,
including local/global cross terms. Keep mesh-refinement and independent
physical checks separate from algebraic completeness checks.

Acceptance counts: solid torus +1 volume mode; spherical cavity +0 volume modes;
closed torus surface +2 surface modes; annulus +0 added surface modes with its
existing hole coefficient retained; genus-g surface +2g added surface modes.
For disconnected meshes, sum component counts and remove gauges per component.
