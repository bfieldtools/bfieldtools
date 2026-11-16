# Geometry, topology, and matrix assembly

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

## Resistance and planned inductance

For component operators $B_k$ with shape $(N_T,N_b)$,

$$
J_k=B_kx,\qquad
R=\sum_{k=x,y,z}B_k^{\mathsf T}\operatorname{diag}(\rho_TV_T)B_k.
$$

Resistance depends on geometry through cell volumes and on the supplied current
operators. It does not need the origin of the coefficient space.

The same representation supports the proposed inductance factorization

$$
M=\frac{\mu_0}{4\pi}\sum_{k=x,y,z}B_k^{\mathsf T}KB_k,
\qquad
K_{TU}=\int_T\int_U\frac{dV_{\mathbf y}\,dV_{\mathbf x}}
{|\mathbf x-\mathbf y|}.
$$

The scalar interaction kernel depends on geometry alone. Its evaluator should
accept mesh geometry and cell-pair indices, without knowing edge ordering,
topology, gauge, resistivity, or basis normalization. A computed interaction can
then be reused for different current bases on the same geometry.

Inductance integration is not implemented by this refactor. When it is added,
the contraction can process cell blocks so that storing a full dense cell-pair
matrix is optional. The division between geometric integration and current
operators applies to dense, blocked, and matrix-free evaluation.

The volume factors belong inside $K$. Since $B_k$ already contains current
density, the inductance contraction must not multiply by cell volumes again.
The kernel derivations and quadrature limitations are documented in
[INDUCTANCE_INTEGRALS.md](INDUCTANCE_INTEGRALS.md).
