# Volumetric eddy-current modeling

Status: proposed plan for review. No implementation decisions below are final.

## Objective and initial scope

Extend bfieldtools from surface currents to divergence-free volume currents on
tetrahedral meshes. Assemble a sparse resistance matrix and a dense inductance
matrix, then solve

$$
R v_k = \lambda_k M v_k
$$

to obtain current modes with decay times $\tau_k = 1/\lambda_k$. Use these modes
for driven and freely decaying eddy-current simulations, similarly to the existing
surface workflow.

The first implementation assumes:

- Isolated conductors with no electrical terminals or imposed transport current.
- Linear, time-independent, finite conductivity, initially scalar and constant per
  tetrahedron.
- Nonmagnetic materials in free space, with permeability $\mu_0$.
- The magnetoquasistatic approximation: neglect displacement current and retardation.
- Conforming, first-order tetrahedra and the lowest-order compatible current space.

Only conducting regions need a volume mesh. Cavities and insulating gaps must be
represented correctly, but do not need to be filled with computational cells.
Magnetic materials, terminals, and higher-order elements are later extensions.

## First milestone: eddy currents driven by a nearby current loop

Deliver one reproducible example that starts from a closed triangular boundary
mesh of a simply connected conductor and computes its response to a nearby
current loop with a prescribed time-dependent current. A coarse solid sphere
with uniform scalar conductivity is the initial geometry.

The example must:

1. Generate and validate a tetrahedral volume mesh from the boundary triangles.
2. Construct conserved volume-current functions with zero boundary flux and
   independent, gauge-fixed edge coefficients.
3. Assemble the sparse resistance matrix $R$ and dense inductance matrix $M$,
   using the analytic source integration described in
   [INDUCTANCE_INTEGRALS.md](INDUCTANCE_INTEGRALS.md).
4. Compute the generalized eigendecomposition $Rv_k=\lambda_kMv_k$, reporting
   current modes and their decay times $1/\lambda_k$.
5. Calculate coupling to an external current loop, use a changing loop current
   to excite the modes, and reconstruct the induced volume-current density.

Place the loop outside and close to the conductor, with no intersection. For
unit loop current, let $\mathbf A_{\mathrm{loop}}^{(1)}$ be its vector potential
and define the coupling vector

$$
m_a=\int_\Omega\mathbf j_a\cdot\mathbf A_{\mathrm{loop}}^{(1)}\,dV.
$$

For prescribed loop current $I(t)$, the excitation and modal equations are

$$
M\dot x+Rx=-m\,\dot I(t),
\qquad
\dot a_k+\lambda_k a_k=-(v_k^{\mathsf T}m)\,\dot I(t).
$$

Use a smooth current pulse with an initially unexcited conductor, followed by
a constant-current or zero-current interval to demonstrate free decay. A static
loop current supplies a magnetic field but no continuing eddy-current excitation
in this stationary-conductor model. The loop current is prescribed; solving the
source circuit and its back reaction is outside this milestone.

Acceptance requires a plot of the input waveform, mode decay times and amplitudes,
and volume-current slices or vectors at selected times. Check conservation,
boundary leakage, matrix symmetry and positive energy, eigenpair residuals, and
convergence of representative low modes and coupling coefficients. Compare the
modal response with direct integration of the same discrete system on the coarse
mesh. Broader geometry/topology validation and performance optimization follow
this first complete example.

## Geometry processing and meshing

Use a backend-independent tetrahedral mesh container. Keep meshing dependencies
optional so that users can also supply an existing tetrahedral mesh.

| Tool | Role | Proposed use |
| --- | --- | --- |
| [TetGen through Python `tetgen`](https://tetgen.pyvista.org/api.html) | Tetrahedralization of closed triangle meshes; holes, material regions, and size constraints | First meshing adapter; accepts vertex and triangle arrays directly |
| [Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#t13) | CAD and discrete-surface geometry, local size control, multiple regions | Alternative backend for broader geometry workflows; tutorial `t13` covers STL remeshing |
| [CGAL / pygalmesh](https://github.com/meshpro/pygalmesh#volume-meshes-from-surface-meshes) | Quality meshing from surfaces, implicit geometry, and labelled volumes | Evaluate if those input types become important |
| [fTetWild](https://github.com/wildmeshing/fTetWild) | Meshing imperfect triangle soups with a geometric tolerance | Fallback for difficult input surfaces; verify preservation of thin gaps and topology |
| [trimesh](https://trimesh.org/trimesh.repair.html) | Surface loading, inspection, and simple repairs | Reuse the existing dependency for input preparation |
| PyVista | Tetrahedral visualization, slicing, and current display | Reuse the existing dependency |

Input validation should check closure, manifold connectivity, consistent
orientation, degenerate triangles, self-intersections, and nested components.
An outer boundary alone is insufficient when a conductor has internal cavities or
material interfaces. Repairs that alter geometry or topology should be explicit.

Meshing must preserve boundary and material labels. Validate positive cell volumes,
face conformity, cell quality, and the intended number of connected components and
holes. Refinement should eventually resolve the skin depth and relevant geometric
features, rather than use a single global cell size exclusively.

## Current representation

### From loops around surface vertices to loops around volume edges

This section builds on the surface-current construction in
[Mäkinen et al., Part I, Sections II–III](https://doi.org/10.1063/5.0016090).
A surface hat function generates a current loop around its vertex. In a
tetrahedral volume mesh, the corresponding local building block is a current
loop around an **interior edge**.

Imagine the tetrahedra incident on an interior edge as wedges around an axis.
The basis current passes from wedge to wedge, circulating around the edge.
Within each tetrahedron its direction and magnitude are constant. Across the
shared faces, the normal current matches, so the wedges carry one continuous
circulating pattern. The edge identifies the axis of this pattern; current is
spread through the surrounding cells rather than carried along the edge itself.

We generate these patterns with a current vector potential $\mathbf T$:

$$
\mathbf T_h=\sum_e t_e\mathbf N_e,
\qquad
\mathbf J_h=\nabla\times\mathbf T_h
=\sum_e t_e\mathbf j_e,
\qquad
\mathbf j_e=\nabla\times\mathbf N_e.
$$

Here $\mathbf N_e$ is the potential function associated with edge $e$, and
$\mathbf j_e$ is its circulating current pattern. The potential $\mathbf T$
parametrizes current and is distinct from the magnetic vector potential
$\mathbf A$ used for field coupling.

| Surface formulation | Volume counterpart |
| --- | --- |
| Scalar stream function $\psi$ | Current vector potential $\mathbf T$ |
| Linear nodal hat function | Whitney/Nédélec edge function |
| Rotated gradient gives current | Curl gives current |
| Constant current within each triangle | Constant current within each tetrahedron |
| Local loop around a vertex | Local loop around an interior edge |

### What current does one edge coefficient represent?

The potential functions are normalized by their line integrals along mesh edges:

$$
\int_{e'}\mathbf N_e\cdot d\mathbf l=\delta_{e'e},
\qquad
 t_e=\int_e\mathbf T_h\cdot d\mathbf l.
$$

Although this is an integral of a potential, its units are amperes:
$\mathbf T$ has units A/m. Its current interpretation follows from taking a
cross section through the loop.

Excite only one interior edge $e=(i,j)$, so that
$\mathbf J_h=t_e\mathbf j_{ij}$. Choose any triangular face $f$ containing that
edge. Such a face cuts across the circulating pattern as it passes between two
neighbouring tetrahedra. Stokes' theorem gives

$$
\int_f\mathbf J_h\cdot\mathbf n_f\,dS
=t_e\oint_{\partial f}\mathbf N_{ij}\cdot d\mathbf l
=\pm t_e.
$$

The other two edges of the face contribute zero to the line integral. Thus
**for this single basis pattern, $t_e$ is the total current circulating around
edge $e$**, measured through a face that cuts the loop. Orientation determines
the sign. Each face containing $e$ measures the same magnitude; these are
successive cross sections of one circulation, not separate parallel currents.

For a superposition of patterns, the measured face current is instead

$$
q_f=\int_f\mathbf J_h\cdot\mathbf n_f\,dS
=\sum_{e\subset\partial f}C_{fe}t_e,
$$

where $C_{fe}$ is the signed face–edge incidence matrix. Several edge patterns
can contribute to the same cross section. This distinction matters when
interpreting individual coefficients of a computed mode.

### Stream surfaces explain the shape of an edge loop

The local geometry is especially clear in barycentric coordinates. Let $\ell_i$
and $\ell_j$ be the linear hat functions at the endpoints of the edge. The
[Whitney/Nédélec edge function](https://defelement.org/elements/nedelec1.html)
and its current are

$$
\mathbf N_{ij}=\ell_i\nabla\ell_j-\ell_j\nabla\ell_i,
\qquad
\mathbf j_{ij}=2\nabla\ell_i\times\nabla\ell_j.
$$

The current is perpendicular to both gradients. Its streamlines therefore lie
on the intersections of the two families of surfaces

$$
\ell_i=\text{constant},\qquad \ell_j=\text{constant}.
$$

Within a tetrahedron these are planes, and their intersection is a straight line
parallel to the edge opposite $ij$. Across the tetrahedra around an interior
edge, the segments join into polygonal loops. The barycentric surfaces provide
two labels for these loops, playing the geometric role that stream-function
contours play on a surface.

This is a direct connection to the **Euler-potential representation**, closely
related to Clebsch variables:

$$
\mathbf J=\nabla\alpha\times\nabla\beta.
$$

Its current lines follow intersections of $\alpha$ and $\beta$ level surfaces.
For an individual edge pattern, the potentials are simply barycentric functions,
with the factor $2t_e$ supplying the amplitude. For an introduction to the
stream-surface and flux interpretation, see
[Euler potentials for two line currents (2011)](https://doi.org/10.1029/2011JA016595).

The normalization also has a geometric interpretation. For the weighted edge
pattern, the oriented current through a small bundle of streamlines is

$$
dI=2t_e\,d\ell_i\,d\ell_j.
$$

On a triangular face containing $ij$, the coordinate region is
$\ell_i\geq0$, $\ell_j\geq0$, $\ell_i+\ell_j\leq1$. Its area in the
$(\ell_i,\ell_j)$ plane is $1/2$, giving total current $\pm t_e$, consistent
with Stokes' theorem.

These stream surfaces describe one basis pattern. A general linear combination
need not share one global pair of stream-surface labels. A single smooth,
single-valued Euler-potential pair restricts field-line topology; see
[Price (2010)](https://academic.oup.com/mnras/article/401/3/1475/1089493).
We retain a linear sum of local patterns for the matrix formulation rather than
require a global two-scalar parametrization.

### Physical fluxes are unique; the decomposition into edge loops is not

As in the surface formulation, a potential representation has freedom that does
not change current. In the volume case,

$$
\mathbf T_h\longmapsto\mathbf T_h+\nabla\psi_h,
\qquad
 t_{ij}\longmapsto t_{ij}+\psi_j-\psi_i
$$

leaves $\mathbf J_h$ and every face flux unchanged. Consequently, an individual
$t_e$ cannot be recovered uniquely from the total current field without choosing
a **gauge**. Its interpretation as a loop amplitude is exact for an isolated
basis pattern, but amplitudes in a redundant superposition can be redistributed
without changing any physical current.

The gauge-independent relation is the closed-contour integral

$$
\oint_{\partial S}\mathbf T_h\cdot d\mathbf l
=\int_S\mathbf J_h\cdot\mathbf n\,dS.
$$

An open-edge integral has this interpretation only after specifying a return
path and accounting for its contribution. This is analogous to stream-function
differences, rather than absolute stream-function values, measuring current on
a surface.

For computation, remove the redundant combinations before solving the generalized
eigenproblem. The volume gradient freedom is larger than the single additive
constant of a scalar surface stream function; deleting one coefficient is
usually insufficient.

The proposed solver unknowns are **gauge-fixed independent edge coefficients**.
For a simply connected conductor, use a boundary-compatible **tree gauge**:
choose a suitable spanning tree and use the nodal potential $\psi_h$ to set its
edge coefficients to zero. On an unconstrained connected mesh, this follows by
solving $\psi_j-\psi_i=-t_{ij}$ successively along the tree. With boundary
conditions, adapt the tree and allowed potential values to preserve the boundary
constraints; do not choose and prune tree edges independently of those constraints.

For a retained edge, its return path through the zero-coefficient tree forms a
closed contour $\gamma_e$. When it bounds a surface $S_e$ in the conductor,
consistent orientation gives

$$
t_e=\oint_{\gamma_e}\mathbf T_h\cdot d\mathbf l
=\int_{S_e}\mathbf J_h\cdot\mathbf n\,dS.
$$

Thus the retained coefficient measures current through a surface spanning that
contour. After enforcing boundary conditions and removing the full zero-current
freedom, these coefficients uniquely parametrize the represented current density.
Assemble $R$ and $M$ directly in this independent edge basis. For conductors with
handles, also address curl-free freedoms beyond scalar gradients and include the
required global current patterns, as discussed below. Face fluxes and dual cycles
serve as reconstruction and validation tools; they need not be solver unknowns.

### Face basis functions reconstruct the same volume currents

The edge-loop picture supplies physical building blocks. Face fluxes provide a
second description that makes reconstruction, conservation, and completeness easy
to check.

Use the **lowest-order Raviart–Thomas space**, $RT_0$, to reconstruct a field from
one oriented flux $q_f$ per face. Its local fields have the form
$\mathbf a_T+b_T\mathbf r$, with four parameters determined by the four face
fluxes. See the [Raviart–Thomas element definition](https://defelement.org/elements/raviart-thomas.html);
some libraries call this degree one.

For a tetrahedron $T$ of volume $V_T$, the local function for outward unit flux
through the face opposite vertex $\mathbf r_i$ is

$$
\boldsymbol\phi_{T,i}(\mathbf r)=\frac{\mathbf r-\mathbf r_i}{3V_T}.
$$

It has zero flux through the other three faces. Join the two local functions at
an interior face with signs consistent with its chosen normal to obtain
$\boldsymbol\phi_f$. Then

$$
\mathbf J_h=\sum_f q_f\boldsymbol\phi_f,
\qquad
\int_g\boldsymbol\phi_f\cdot\mathbf n_g\,dS=\delta_{gf}.
$$

Thus a face coefficient is the current crossing that face, while the associated
basis function specifies how that current is distributed in its two neighbouring
cells. A single face function has an imbalance in those cells. Conservation
selects combinations whose fluxes balance.

Let $D$ be the signed cell–face incidence matrix, with signs given by outward
versus inward face orientation. For an RT0 field,

$$
(\nabla\cdot\mathbf J_h)|_T=\frac{(Dq)_T}{V_T}.
$$

Setting $Dq=0$ forces $b_T=0$, since
$\nabla\cdot(\mathbf a_T+b_T\mathbf r)=3b_T$. The conserved fields are therefore
constant in each cell, with matching normal components, exactly as in the
edge-loop construction. Stokes' theorem connects their coefficients by $q=Ct$;
$DC=0$ expresses that each edge loop conserves current.

### The function spaces behind the geometry

The face space is **$H(\mathrm{div})$ conforming**: the field and its weak
divergence are square-integrable. For these piecewise-polynomial fields, this
means the normal component matches across faces, while tangential components may
jump. It prevents artificial interface sources but does not alone impose zero
divergence. See [Schöberl's introduction to $H(\mathrm{div})$](https://jschoeberl.github.io/iFEM/secondorder/hdiv.html).

The potential edge space is **$H(\mathrm{curl})$ conforming**: the field and its
weak curl are square-integrable. Its tangential components match across faces,
so taking a curl introduces no artificial sheet contributions. Its curls belong
to the RT0 current space. The compatible spaces form the **discrete de Rham
sequence**:

$$
P_1\xrightarrow{\nabla}\mathrm{Nedelec}_0
\xrightarrow{\nabla\times}RT_0
\xrightarrow{\nabla\cdot}P_0.
$$

Here $P_1$ is the continuous linear nodal space and $P_0$ the discontinuous
cellwise-constant scalar space. The sequence places nodal potential differences,
edge integrals, face fluxes, and cell balances in one construction. See
[Schöberl's introduction to the de Rham complex](https://jschoeberl.github.io/iFEM/Maxwell/deRham.html)
for the function-space interpretation. Topology and boundary conditions determine
whether every field with zero derivative comes from the preceding space.

### Boundaries, global loops, and the implementation basis

For an isolated conductor, require zero flux through every boundary face,
including cavity walls. Setting boundary-edge potential coefficients to zero is
one sufficient way to enforce this. However, local loops around interior edges
may then miss currents circulating around a handle of the conductor. Those need
additional global patterns. Conservation of every basis function does not by
itself establish completeness.

A complete reference construction starts directly with face fluxes. Set boundary
fluxes to zero and call the remaining conservation matrix $D_i$. Choose independent
columns of $Z$ spanning $\ker D_i$:

$$
q_i=Zx,\qquad
\mathbf j_a=\sum_{f\ \mathrm{interior}}Z_{fa}\boldsymbol\phi_f.
$$

These $\mathbf j_a$ form an independent conserved-current basis for $R$ and $M$.
One way to construct $Z$ is to draw the **dual graph**: tetrahedra are nodes and
interior faces connect neighbouring nodes. Each independent graph cycle routes
current through a closed sequence of cells. A spanning forest gives one cycle
per remaining graph edge, with count

$$
N_{\mathrm{cycles}}=N_{\mathrm{interior\ faces}}-N_{\mathrm{tetrahedra}}
+N_{\mathrm{dual\ components}}.
$$

This includes circulation around handles automatically, although long cycles
increase basis overlap and fill in $R$.
[Passarotto, Pitassi and Specogna](https://www.comphys.com/pdf/cmame_22_vinco.pdf)
compare loop and vector-potential constructions and their treatment of topology.

Use face reconstruction and dual cycles as the reference for validating the
proposed gauge-fixed edge basis, supplemented with global patterns where required.
Both should represent the same admissible currents. Check locality and conditioning
before committing to a particular tree construction.

## Resistance, inductance, and coupling

For independent current functions $\mathbf j_a$, define

$$
R_{ab} = \int_\Omega
\mathbf j_a(\mathbf r)\cdot\boldsymbol\rho(\mathbf r)\mathbf j_b(\mathbf r)\,dV,
$$

$$
M_{ab} = \frac{\mu_0}{4\pi}
\int_\Omega\int_\Omega
\frac{\mathbf j_a(\mathbf r)\cdot\mathbf j_b(\mathbf r')}
{|\mathbf r-\mathbf r'|}\,dV'\,dV.
$$

Here $\boldsymbol\rho$ is resistivity; the initial implementation uses scalar
$\rho=1/\sigma$. The physical quantities are Joule loss
$x^\mathsf TRx$ and magnetic energy $\tfrac12 x^\mathsf TMx$.

In current-vector-potential coefficients, resistance is a resistivity-weighted
curl–curl matrix rather than a scalar Laplacian. Local current support makes
resistance sparse; inductance is nonlocal and dense. Both matrices must use the
same independent current space.

For lowest-order conserved RT currents, the current is constant within each
tetrahedron. Let $P_\alpha$ map coefficients to the cellwise current-density
component $\alpha$, and define

$$
K_{tu} = \int_{T_t}\int_{T_u}
\frac{1}{|\mathbf r-\mathbf r'|}\,dV'\,dV.
$$

For cellwise scalar resistivity, this gives

$$
R = \sum_{\alpha=x,y,z} P_\alpha^\mathsf T
\operatorname{diag}(\rho_t V_t)P_\alpha,
\qquad
M = \frac{\mu_0}{4\pi}\sum_{\alpha=x,y,z}
P_\alpha^\mathsf T K P_\alpha.
$$

This separates current reconstruction from the scalar tetrahedron interaction
kernel and resembles the existing surface assembly architecture. Whether to store
the cell interaction matrix or assemble reduced inductance directly should depend
on memory estimates for the selected mesh and basis.

Treat far, near, touching, and identical tetrahedron pairs separately. Standard
quadrature is suitable for well-separated cells; self and touching interactions
need singularity-aware integration. Do not substitute arbitrary diagonal values or
evaluate the singular kernel at coincident quadrature points.

Add vector-potential and magnetic-field coupling operators alongside inductance.
They support both excitation and inspection of computed current modes.

## Modes and time evolution

After removing redundant coefficients, solve

$$
Rv_k=\lambda_kMv_k, \qquad v_i^\mathsf TMv_j=\delta_{ij}.
$$

For the initial passive model, eigenvalues are positive decay rates, with
$\tau_k=1/\lambda_k$. Begin with a dense symmetric generalized eigensolver on
small problems.

For an externally imposed vector potential,

$$
M\dot x+Rx=f(t), \qquad
f_a(t)=-\int_\Omega\mathbf j_a\cdot\partial_t\mathbf A_{\mathrm{ext}}\,dV.
$$

With $x=\sum_k v_k a_k$, modal evolution is

$$
\dot a_k+\lambda_k a_k=v_k^\mathsf Tf(t).
$$

Check reduced simulations against direct integration of the full discrete system.
For free decay, verify the energy identity

$$
\frac{d}{dt}\left(\tfrac12x^\mathsf TMx\right)=-x^\mathsf TRx.
$$

## Implementation sequence and acceptance criteria

Steps 1–8 deliver the first milestone for the simple conductor and nearby loop.
Their broader geometry and topology examples are follow-up coverage, rather than
prerequisites for that first demonstration.

1. **Mesh representation and geometry contract.** Store vertices, oriented
   tetrahedra, faces and edges, adjacency, boundary tags, material labels, cell
   volumes, and barycentric gradients. Accept precomputed meshes. Verify geometry
   and orientation on small hand-constructed examples.

2. **First meshing adapter.** Implement closed triangle mesh → TetGen → validated
   volume mesh. Provide a sphere, block, and solid torus example, plus a cavity
   example. Inspect cell quality and boundary preservation with PyVista slices.

3. **Topology and current reconstruction.** Assemble incidence operators, eliminate
   boundary fluxes, and construct the independent dual-cycle reference basis.
   Implement local edge curls and a boundary-compatible tree gauge, retaining
   independent edge coefficients as the proposed solver unknowns. Account for
   additional topology-dependent freedoms and global current patterns.
   Check conservation, zero
   boundary leakage, rank, expected cycle count, and current-space equivalence.
   Include disconnected conductors and toroidal circulation.

4. **Resistance.** Assemble for scalar resistivity constant per cell. Check
   symmetry, positive definiteness in independent coefficients, and agreement
   with directly integrated Joule loss. Visualize reconstructed cell currents.
   This is an intermediate checkpoint toward the first milestone.

5. **Tetrahedron interaction kernel.** Establish a convergent reference integration
   method for self, touching, near, and far pairs. Check symmetry, geometric
   scaling, and convergence before optimizing assembly.

6. **Inductance and field coupling.** Assemble small dense problems. Verify positive
   magnetic energy, reciprocity, and quadrature convergence. Implement vector
   potential and magnetic-field evaluation, documenting accuracy near and inside
   cells. Compare reconstructed-current evaluations with independent integration.

7. **Eigenmodes.** Solve the generalized problem and inspect current patterns and
   decay times. Check residuals and M-orthogonality. Confirm invariance under basis
   changes and convergence of low modes under mesh and quadrature refinement.
   Compare invariant subspaces when eigenvalues are degenerate.

8. **Excitation and transients.** Calculate the nearby-loop coupling vector and
   modal excitation from a prescribed current pulse. Reconstruct induced volume
   currents during excitation and free decay. Compare modal evolution with
   full-system time integration, verify energy balance, and measure
   modal-truncation error. Complete the reproducible first-milestone example.

9. **Physical validation.** Compare suitable sphere/slab responses with analytical
   solutions and a thin conductor with the existing surface formulation. Exercise
   skin-depth refinement, multiple components, cavities, and handles. Distinguish
   geometry, discretization, quadrature, and modal-truncation errors.

10. **Performance and broader integration.** Profile assembly and eigensolves;
    introduce compression or matrix-free inductance only after the reference
    implementation is validated. Add higher-level conductor APIs and examples
    once the low-level representation is stable.

## Proposed module boundaries

Names are provisional; keep the first implementation small.

| Module | Responsibility |
| --- | --- |
| `mesh.py` | Tetrahedral mesh container, geometry, labels, validation |
| `meshing.py` | Optional meshing adapters |
| `topology.py` | Incidence matrices, boundary constraints, cycle bases and gauges |
| `basis.py` | RT face functions, edge curls, independent current reconstruction |
| `impedance.py` | Resistance and inductance assembly |
| `integrals.py` | Tetrahedron interaction kernels and singular integration |
| `magnetics.py` | Vector-potential and magnetic-field coupling |
| `modes.py` | Generalized modes and modal dynamics |

Reuse existing numerical and visualization conventions where practical, while
keeping volume-current operators distinct from surface stream-function operators.
Packaging currently explicitly lists only `bfieldtools`; update package inclusion
when the new `bfieldtools.volume` subpackage becomes usable.

## Literature and reference implementations

- **Passarotto, Pitassi and Specogna (2022),**
  [Foundations of volume integral methods for eddy current problems](https://www.comphys.com/pdf/cmame_22_vinco.pdf).
  Main implementation reference: conserved-current formulations, topology,
  uniform-volume representations, and inductance factorization.
- **Albanese and Rubinacci (1988),**
  [Integral formulation for 3D eddy-current computation using edge elements](https://doi.org/10.1049/ip-a-1.1988.0072).
  Classic current-vector-potential and edge-element formulation.
- **Li and Chew (2006),**
  [Applying divergence-free condition in solving the volume integral equation](https://doi.org/10.2528/PIER05061303).
  Independent and redundant volume-loop constructions and conditioning tradeoffs.
- **Arnold, Falk and Winther (2009/2010),**
  [Finite element exterior calculus: from Hodge theory to numerical stability](https://arxiv.org/abs/0906.4325).
  Mathematical foundation for compatible spaces, kernels, and topology.
- **DefElement:** [Nédélec first kind](https://defelement.org/elements/nedelec1.html)
  and [Raviart–Thomas](https://defelement.org/elements/raviart-thomas.html).
  Reference definitions and degree-of-freedom conventions.
- **[scikit-fem](https://scikit-fem.readthedocs.io/en/stable/api.html):** tetrahedral
  Nédélec and Raviart–Thomas elements as an independent assembly reference;
  not necessarily a core dependency.

## Decisions to review

- TetGen as the first optional meshing backend.
- Gauge-fixed independent edge coefficients as solver unknowns, with a
  boundary-compatible tree gauge and global patterns where needed; face-flux
  topology and dual cycles provide reconstruction and a completeness reference.
- Piecewise-constant conserved currents as the first approximation order.
- Nonmagnetic isolated conductors as the initial physical scope.
- First milestone: boundary-to-volume meshing, $R$ and $M$ assembly, generalized
  eigendecomposition, and eddy-current excitation by a nearby driven current loop.
- Completing geometry, conservation, current visualization, and resistance before
  implementing singular inductance integration.
