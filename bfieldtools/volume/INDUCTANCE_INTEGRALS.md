# Analytic reductions for volume inductance integrals

Status: analytic source reduction and a converged outer-quadrature implementation
are available in [integrals.py](integrals.py); further closed-form reductions
remain investigational.
This note complements [PLAN.md](PLAN.md) and assumes familiarity with the
surface-current formulation in [bfieldtools Part I](https://doi.org/10.1063/5.0016090).

## The geometric integral we need

For the lowest-order conserved volume-current representation, each current basis
function is constant within each tetrahedron. Consequently, the inductance entry
can be assembled as

$$
M_{ab}=\frac{\mu_0}{4\pi}\sum_{T,U}
\left(\mathbf j_a|_T\cdot\mathbf j_b|_U\right)K_{TU},
$$

where the geometric interaction is

$$
K_{TU}=\int_T\int_U\frac{1}{|\mathbf x-\mathbf y|}
\,dV_{\mathbf y}\,dV_{\mathbf x}.
$$

Thus we can calculate the scalar interaction once per tetrahedron pair and reuse
it for all current components. This relies on working with conserved currents;
unconstrained RT0 face functions also contain affine terms.

There are two distinct analytic-integration questions:

1. Can we integrate over a source tetrahedron exactly at an arbitrary observation
   point? **Yes**, and we can reduce it to the existing triangle potential.
2. Can we evaluate the complete interaction of two arbitrary tetrahedra with a
   finite closed-form expression? The reductions below do not establish that.
   Analytic inner integration plus numerical outer integration is already a
   practical approach, analogous to the existing surface implementation.

## Exact source integration using triangle potentials

Define the potential of a uniformly filled source tetrahedron $U$:

$$
u_U(\mathbf x)=\int_U\frac{1}{|\mathbf y-\mathbf x|}\,dV_{\mathbf y}.
$$

Write $R=|\mathbf y-\mathbf x|$. Away from coincidence,

$$
\nabla_{\mathbf y}\cdot\left(\frac{\mathbf y-\mathbf x}{R}\right)
=\frac{3}{R}-\frac{|\mathbf y-\mathbf x|^2}{R^3}
=\frac{2}{R}.
$$

The divergence theorem therefore gives

$$
u_U(\mathbf x)=\frac12\int_{\partial U}
\frac{(\mathbf y-\mathbf x)\cdot\mathbf n}{R}\,dS_{\mathbf y}.
$$

For each planar face $f$, choose any point $\mathbf y_f$ on its plane and let
$\mathbf n_f$ be its outward unit normal. The factor

$$
d_f(\mathbf x)=(\mathbf y_f-\mathbf x)\cdot\mathbf n_f
$$

is the signed distance from $\mathbf x$ to the face plane, positive on the
interior side for the outward-normal convention used here. It is constant over
the face. Define the uniform triangle potential

$$
S_f(\mathbf x)=\int_f\frac{1}{|\mathbf y-\mathbf x|}\,dS_{\mathbf y}.
$$

The exact volume potential is then

$$
\boxed{
u_U(\mathbf x)=\frac12\sum_{f\subset\partial U}
d_f(\mathbf x)S_f(\mathbf x).
}
$$

A tetrahedron needs just four triangle-potential evaluations. The triangle
integrals are already implemented in
[`triangle_potential_uniform`](../integrals.py), using edge logarithms and solid
angles. With its displacement convention $\mathbf x-\mathbf y$, take care to
implement the signed height as defined above rather than infer its sign from
an existing distance helper.

### Interior points and boundary limits

This reduction also holds when $\mathbf x$ is inside $U$. To justify applying
the divergence theorem there, remove a ball of radius $\epsilon$ around
$\mathbf x$. The vector field $(\mathbf y-\mathbf x)/R$ has unit magnitude,
so its flux through the small sphere is of order $\epsilon^2$ and vanishes.
There is no additional point-source contribution.

The potential is finite and continuous on faces, edges, and vertices, so the
formula extends there by limits. This mathematical statement does not guarantee
that an existing floating-point triangle routine evaluates every limiting case
accurately. In particular, terms with $d_f=0$ require care if intermediate edge
logarithms diverge, to avoid expressions such as zero times infinity.

For an interior point, all $d_f$ are nonnegative. Outside the tetrahedron, signed
face contributions can cancel. At large distances this cancellation can lose
precision even though the formula is exact.

### Relation to the literature

Analytic polyhedral point potentials are established results, rather than a new
physical approximation. Relevant sources include:

- Fabbri (2008), [Magnetic flux density and vector potential of uniform polyhedral
  sources](https://doi.org/10.1109/TMAG.2007.908698): elementary expressions for
  uniform polyhedral sources. The author's university provides an
  [abstract and bibliographic record](https://cris.unibo.it/handle/11585/49427).
- Carley and Angioni (2017), [Exact integration of surface and volume
  potentials](https://purehost.bath.ac.uk/ws/files/146284153/jem15_revision_ii.pdf):
  analytic potential integration over elements.
- Wilton et al. (1984), [Potential integrals for uniform and linear source
  distributions on polygonal and polyhedral
  domains](https://doi.org/10.1109/TAP.1984.1143304): foundational potential-integral
  literature.

The short derivation above makes the link to the existing bfieldtools triangle
kernel explicit and can serve as an in-context appendix derivation.

## Outer tetrahedron integration

Once the source potential is known analytically,

$$
K_{TU}=\int_T u_U(\mathbf x)\,dV_{\mathbf x}
\approx V_T\sum_q w_q u_U(\mathbf x_q),
\qquad \sum_q w_q=1.
$$

This replaces a six-dimensional quadrature of the singular kernel with a
three-dimensional quadrature of a finite potential. It follows the existing
surface assembly pattern: analytic source integration and numerical testing
integration.

The outer integrand is not infinitely smooth across the source boundary.
Self interactions, touching cells, and close separated cells therefore still
need convergence checks and potentially subdivision or adaptive quadrature.
Eliminating the explicit singularity does not make a fixed low-order outer rule
exact.

[Passarotto, Pitassi and Specogna (2022), Section 6.2](https://www.comphys.com/pdf/cmame_22_vinco.pdf)
uses analytic inner integration for uniform-volume interactions and compares it
with direct double quadrature. It also constructs a fully analytic reference
for the self interaction of a particular tetrahedron using divergence-theorem
reductions and symbolic computation. That special-case reference does not by
itself provide a general closed-form tetrahedron-pair routine.

Exact interactions satisfy $K_{TU}=K_{UT}$. Numerical outer quadrature may break
this equality. Evaluating both directions is a useful accuracy check; averaging
them can restore symmetry, but cannot substitute for convergence or guarantee
positive definiteness of the assembled matrix.

## Exact reduction of the double volume integral to face pairs

A second derivation removes both volume integrations. For
$R=|\mathbf x-\mathbf y|$,

$$
\Delta_{\mathbf x}R=\frac{2}{R},
\qquad
\nabla_{\mathbf x}R=-\nabla_{\mathbf y}R.
$$

Apply the divergence theorem first in $T$:

$$
K_{TU}=\frac12\int_{\partial T}\int_U
\mathbf n_T\cdot\nabla_{\mathbf x}R\,dV_{\mathbf y}\,dS_{\mathbf x}.
$$

Replace the gradient by $-\nabla_{\mathbf y}R$ and apply the theorem in $U$:

$$
K_{TU}=-\frac12\int_{\partial T}\int_{\partial U}
(\mathbf n_T\cdot\mathbf n_U)R\,dS_{\mathbf y}\,dS_{\mathbf x}.
$$

Since each face normal is constant,

$$
\boxed{
K_{TU}=-\frac12
\sum_{f\subset\partial T}\sum_{g\subset\partial U}
(\mathbf n_f\cdot\mathbf n_g)
\int_f\int_g|\mathbf x-\mathbf y|\,dS_{\mathbf y}\,dS_{\mathbf x}.
}
$$

All normals here are outward normals of the individual tetrahedra. A shared
face consequently has opposite normals for the two cells.

The six-dimensional integral becomes sixteen four-dimensional triangle-pair
integrals. The kernel is now distance $R$, which is finite at coincidence, though
its derivatives are not smooth there. The identity also holds for self and
touching interactions: $\Delta R=2/R$ has no additional delta contribution at
coincidence, and the small excluded-boundary terms vanish.

This is an exact reduction, not yet a fully integrated formula. It suggests a
second numerical evaluator or further analytic reduction using triangle distance
integrals. Cancellation between signed face-pair contributions can again be
significant for widely separated cells.

For related double-divergence reductions of electromagnetic volume interactions,
see [Rivero et al. (2019), 6-D Integrals for Numerical Evaluation by Double
Application of the Divergence Theorem](https://usncursi.org/archive/aps-ursi/2019/abstracts/1962.pdf).
The specific scalar formula above follows from the two identities shown here.

## What is established and what remains open

- The tetrahedron point potential can be evaluated analytically through four
  existing triangle potentials.
- The full pair interaction can be reduced exactly to face-pair integrals of
  distance.
- A finite closed-form evaluator for every arbitrary tetrahedron pair has not
  been established in this investigation. Further recurrence relations or
  special-case derivations would need to be developed or verified from sources.
- Some literature on analytic mutual polyhedral potentials uses infinite
  expansions. For example, [Werner and Scheeres (2005), Mutual Potential of
  Homogeneous Polyhedra](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/42570/10569_2004_Article_4621.pdf)
  explicitly presents a series; it should not be treated as a general finite
  closed-form formula for touching tetrahedra.

## Implemented numerical evaluator

`tetrahedron_potential(mesh, points, source_cells=None)` evaluates source point
potentials using the four-face formula above. Source coordinates are translated
and normalized before using the existing triangle kernel, controlling its fixed
edge-logarithm regularization. Interior, face, edge, and vertex limits are tested.
Beyond eight source radii, a regular 64-point source rule avoids cancellation in
the signed face sum; this is numerical integration of a smooth kernel.

`tetrahedron_interaction_block(mesh, rows, columns)` integrates the resulting
finite potentials over testing cells. It uses an interior positive Gauss–Jacobi
product rule on collapsed tetrahedral coordinates. `quadrature_order=n` means
$n^3$ points, with polynomial exactness through total degree $2n-1$.

Outer integration starts at order 4 by default, compares with order 6, then
increases by two for pairs whose successive values differ by more than
`rtol=1e-4`. The default `max_order=12` limits refinement. Nonconvergence raises
an error identifying the pair; increase `max_order` or relax `rtol` explicitly.
The estimate is empirical, not a rigorous error bound. The far-source rule and
floating-point triangle routine also limit the accuracy attainable by tightening
outer quadrature alone.

Both source/testing directions are evaluated and averaged to enforce
reciprocity. `iter_tetrahedron_interactions` yields upper-triangular cell blocks;
`tetrahedron_interaction_matrix` materializes the complete symmetric $K$ when
needed. `inductance_matrix` normally contracts blocks directly into the current
basis, retaining no full cell-pair matrix. It can also reuse a supplied $K$.
Neither evaluator changes diagonal values or eigenvalues to force positivity.

Tests include:

- Positive quadrature weights and exact barycentric moments.
- Point potentials against independent numerical face integration.
- Self interaction against the face-pair distance reduction.
- The summed interaction of a tetrahedral unit cube against a smooth,
  independent displacement-space integral.
- Shared-face, shared-edge, shared-vertex, near, and far cell pairs.
- Reciprocity, rigid motions, fifth-power scaling, and order refinement.
- Blocked versus precomputed-kernel assembly, magnetic energy, basis changes,
  generalized eigenpair residuals, M-orthogonality, and inverse-square scaling
  of decay rates with conductor size.

`examples/volume/eddy_modes.py` shows the lowest generalized modes and provides
`--check-convergence` to repeat integration with a tighter tolerance. Its sphere
mesh is deliberately coarse; quadrature convergence does not establish spatial
mesh convergence or complete physical validation.

## Further integration work

1. Extend conditioning and boundary-limit validation to very thin tetrahedra and
   extreme coordinate scales.
2. Add spatial subdivision or adaptive cell integration when order refinement
   is inefficient or fails on difficult self, touching, and close pairs.
3. Expand high-accuracy independent references beyond the cube and current
   face-pair comparisons, including verified analytic special cases.
4. Optimize far-pair evaluation with controlled quadrature or moment expansions,
   and expose its accuracy independently of outer integration if tighter
   tolerances are needed.
5. Investigate further analytic reductions if profiling shows close-pair assembly
   is a significant cost.

Validation should cover:

- Reciprocity: $K_{TU}=K_{UT}$.
- Translation and rotation invariance.
- Uniform scaling: $K_{sT,sU}=s^5K_{TU}$ for $s>0$.
- Far-field behaviour: $K_{TU}\sim V_TV_U/d$, with centroid separation $d$.
- Self, shared-face, shared-edge, shared-vertex, and close separated pairs.
- Convergence under outer quadrature refinement and cell subdivision.
- Agreement between independent formulations and verified special cases.
- Positive magnetic energy after assembly in the independent current basis.

The implemented baseline is therefore **analytic source integration plus
controlled outer quadrature**, with regular far-source quadrature. The face-pair
identity remains an independent validation route and a starting point for a more
complete analytic derivation.
