# Geometry, topology, and matrix assembly

See [PLAN.md](PLAN.md) for status and remaining work,
[THEORY.md](THEORY.md) for the equivalent-magnetization interpretation and
formulation, and [VALIDATION.md](VALIDATION.md) for tests and example workflows.

The volume implementation uses composition, with four separate data containers.
Algorithms live in functions; constructor helpers and cached incidence properties
delegate to those functions.

| Container | Stores | Responsibility |
| --- | --- | --- |
| `TetraMesh` | Coordinates, oriented cells, volumes, boundary triangles, cached geometry | Geometry and cell ordering |
| `TetraTopology` | Ordered cell indices, vertex count, edge/face connectivity, signs, boundary masks, cached incidence matrices | Geometry-independent discrete topology |
| `TreeCotreeGauge` | Tree and retained edges, coefficient prolongation, reduced face-flux map | A particular independent coefficient selection |
| `CellCurrentBasis` | Mesh and three real CSR component operators | Piecewise-constant Cartesian current representation for numerical assembly |

`TetraTopology` retains no vertex coordinates or mesh reference. Its cell-index
array shares the mesh's read-only connectivity. It can be reused for a translated,
rotated, or scaled geometry with identical ordered connectivity. Gauge selection
likewise depends only on topology.

`CellCurrentBasis` retains no topology or gauge. Its mesh reference establishes
the exact geometry and cell ordering of its operators, without copying mesh data.
Its coefficient count is derived from operator shape. Matrix assembly and
reconstruction never traverse a topology/gauge object chain.

## Construction

```python
topology = TetraTopology.from_mesh(mesh)
gauge = TreeCotreeGauge.from_topology(topology)
basis = build_edge_current_basis(mesh, topology, gauge)

R = resistance_matrix(basis, resistivity)
J = reconstruct_currents(basis, coefficients)
```

The edge builder is the meeting point of geometry, topology, and gauge. It
validates that geometry and topology have identical cell ordering, then builds
component operators from the selected edges. It does not select a gauge silently.
The lower-level `cell_current_operators(mesh, topology, retained_edges)` accepts
an ordered edge selection without knowing how it was chosen.

This replaces the earlier `EdgeCurrentBasis.from_topology(topology)` API. The
edge construction now returns a general `CellCurrentBasis`; callers keep topology
and gauge separately when they need face fluxes or edge-coefficient interpretation.
Operators are computed once during construction and retained directly, rather
than retaining all construction inputs for a lazy computation.

Alternative constructors can provide operators directly:

```python
basis = CellCurrentBasis(mesh, (Bx, By, Bz))
```

This supports changes of basis, future global current patterns, and alternative
conserved-current constructions without changing matrix assembly. The container
checks that the three matrices are real, finite CSR matrices with the same shape
and one row per mesh cell. Independence, conservation, and boundary conditions
must be established by the construction algorithm, not inferred from matrix shape.

## Resistance and inductance

For component operators $B_k$ with shape $(N_T,N_b)$,

$$
J_k=B_kx,\qquad
R=\sum_{k=x,y,z}B_k^{\mathsf T}\operatorname{diag}(\rho_TV_T)B_k.
$$

Resistance depends on geometry through cell volumes and on the supplied current
operators. It does not need the origin of the coefficient space.

The same representation supports the implemented inductance factorization

$$
M=\frac{\mu_0}{4\pi}\sum_{k=x,y,z}B_k^{\mathsf T}KB_k,
\qquad
K_{TU}=\int_T\int_U\frac{dV_{\mathbf y}\,dV_{\mathbf x}}
{|\mathbf x-\mathbf y|}.
$$

The scalar interaction kernel depends on geometry alone. Its evaluator accepts
mesh geometry and cell-pair indices, without knowing edge ordering,
topology, gauge, resistivity, or basis normalization. A computed interaction can
then be reused for different current bases on the same geometry.

`inductance_matrix(basis)` processes cell blocks without storing a full dense
cell-pair matrix. Within each block, dense contraction is restricted to active
basis coefficients. The final inductance matrix is dense. Alternatively,
`inductance_matrix(basis, interactions=K)` reuses a precomputed geometric matrix.
`tetrahedron_interaction_matrix(mesh)` constructs $K$; the lower-level block
iterator supports streaming. Matrix-free evaluation is a later extension.

The volume factors belong inside $K$. Since $B_k$ already contains current
density, the inductance contraction must not multiply by cell volumes again.
The kernel derivations and quadrature limitations are documented in
[INDUCTANCE_INTEGRALS.md](INDUCTANCE_INTEGRALS.md).

## Excitation and evolution

`loop_coupling(basis, loop_vertices)` returns the coupling vector in the basis
coefficient space. It needs only current operators and geometry. No topology or
gauge reference is retained by the coupling or time-integration code.

For M-orthonormal generalized eigenvectors V and positive rates λ, use
`modal_response(rates, V.T @ coupling, current_derivative, times)`. This solves
`da/dt + λ*a = -(V.T @ coupling)*dI/dt`, returning amplitudes at the requested
times. Reconstruct coefficients with `x = V @ a` and cell currents using
`reconstruct_currents(basis, x)`. Optional initial amplitudes specify a nonzero
starting state; time evolution otherwise starts from an unexcited conductor.
