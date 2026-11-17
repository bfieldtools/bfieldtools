# Finite-permeability magnetostatic BEM

Status: stage 1 implemented. `MagneticBoundary` validates complete metal
boundaries and material labels; `single_layer_potential`,
`single_layer_gradient`, and `adjoint_double_layer` provide analytic
constant-face operators, reusing the existing triangle kernels.

Run `python examples/validation/finite_permeability_bem.py` for a solid sphere
and spherical shell in a uniform field, with analytic reference curves. The
example assembles and solves the dense equation directly; a reusable solver
API remains part of stage 2. Boundary gradients currently support off-surface
targets only; explicit one-sided trace evaluation is still to be implemented.

Validation covers independent triangle quadrature, potential gradients,
normal-derivative jumps, scaling, chunking, geometry/material validation,
sphere and shell convergence, and close inner/outer surface interactions.
The wheel includes the new subpackage and has been checked for importability.
At $\mu_r=10$ with 1280 faces per spherical surface, the example gives center
fields of approximately $0.2854\,\mathrm{A/m}$ for the solid sphere and
$0.4763\,\mathrm{A/m}$ for the shell, compared with analytic values of
$0.2500$ and $0.4147\,\mathrm{A/m}$ respectively. These finite-resolution
errors remain despite linear-system residuals near machine precision. Mesh
refinement is required, especially for high permeability and thin shields.

## Objective and scope

Extend bfieldtools shielding calculations from the infinite-permeability
approximation to a finite, uniform permeability in each metal body. Use a
boundary-only discretization of the actual metal volume, including both sides
of a hollow shield and surfaces around openings.

The initial implementation assumes static fields, linear isotropic materials,
positive scalar permeability, vacuum outside the metal, and prescribed free
currents that do not intersect the metal. Different disjoint metal bodies may
have different permeabilities. Nonlinear saturation, hysteresis, anisotropy,
eddy currents, and a zero-thickness shell approximation are outside the initial
scope. Uniform permeability does not imply uniform magnetization.

## 1. Field equations and conventions

Let $\Omega_m$ denote a metal volume with permeability $\mu_m$, and let
$\Gamma=\partial\Omega_m$. The surrounding vacuum has permeability $\mu_0$.
The unit normal $\hat{\mathbf n}$ points out of the metal on every boundary
component. On an inner shield surface, it points into the cavity.

Follow Mäkinen et al., *Magnetic-field modeling with surface currents. Part I*,
Sec. II B and Eq. (4), which define the vacuum scalar potential by
$\mathbf B=-\mu_0\nabla U$. Extend this convention to each homogeneous material:

$$
\boxed{\mathbf H=-\nabla U,\qquad \mathbf B=-\mu\nabla U.}
$$

Here $U$ has units of amperes. The potential is an $\mathbf H$ potential;
the permeability factor belongs in the relation to $\mathbf B$.

In a current-free homogeneous region,

$$
\nabla\times\mathbf H=0,\qquad
\nabla\cdot\mathbf B=0,\qquad
\nabla^2 U=0.
$$

With no free surface current on $\Gamma$, the transmission conditions are

$$
U_m=U_v,\qquad
\mu_m\partial_n U_m=\mu_0\partial_n U_v,
\qquad \partial_n=\hat{\mathbf n}\cdot\nabla.
$$

Both normal derivatives use the same $\hat{\mathbf n}$, rather than separate
outward normals for the two regions. Potential continuity enforces tangential
$\mathbf H$ continuity; the second condition enforces normal $\mathbf B$
continuity. Additive potential constants
are chosen consistently across the interface.

A globally single-valued total scalar potential need not exist around linked
free currents. Implement the more general reduced-potential formulation:

$$
\mathbf H=\mathbf H_0-\nabla u,\qquad
\mathbf B=\mu(\mathbf x)\bigl(\mathbf H_0-\nabla u\bigr),
$$

where $\mathbf H_0$ is the field of the prescribed currents with all metal
replaced by vacuum. It is smooth across $\Gamma$, curl-free near and inside the metal, and
divergence-free. The induced correction carries no additional free-current
circulation. Where a primary potential $U_0$ is available,

$$
\mathbf H_0=-\nabla U_0,\qquad U=U_0+u.
$$

## 2. Single-layer representation

Use the free-space Laplace Green function

$$
G(\mathbf x,\mathbf y)=\frac{1}{4\pi|\mathbf x-\mathbf y|}
$$

and represent the induced potential by a surface density $q$:

$$
u(\mathbf x)=(Sq)(\mathbf x)
=\int_\Gamma G(\mathbf x,\mathbf y)q(\mathbf y)\,\mathrm dS_{\mathbf y}.
$$

The density $q$ has units $\mathrm{A/m}$; $u$ has units $\mathrm A$. The
single-layer potential is continuous across $\Gamma$ and decays at infinity.
It is harmonic off $\Gamma$, so
only the normal-flux transmission condition remains to be enforced.

Define the average normal-derivative operator

$$
(K'q)(\mathbf x)
=\operatorname{p.v.}\!\int_\Gamma
\partial_{n_{\mathbf x}}G(\mathbf x,\mathbf y)q(\mathbf y)\,\mathrm dS_{\mathbf y},
$$

$$
\partial_{n_{\mathbf x}}G(\mathbf x,\mathbf y)
=-\frac{\hat{\mathbf n}(\mathbf x)\cdot(\mathbf x-\mathbf y)}
{4\pi|\mathbf x-\mathbf y|^3}.
$$

The derivative is with respect to the target coordinate $\mathbf x$. It must not be
confused with the source-normal derivative in the double-layer operator.

At a smooth point of $\Gamma$, the single-layer normal derivatives satisfy

$$
\left.\partial_n Sq\right|_m=K'q+\frac12q,\qquad
\left.\partial_n Sq\right|_v=K'q-\frac12q.
$$

These geometric jump relations are independent of the field-potential sign
convention. Because $\mathbf H=\mathbf H_0-\nabla Sq$, they imply

$$
q=\hat{\mathbf n}\cdot(\mathbf H_v-\mathbf H_m)
=\mathbf M\cdot\hat{\mathbf n},\qquad
\mathbf M=(\mu_r-1)\mathbf H_m,
$$

where the second equality assumes vacuum outside and
$\mu_r=\mu_m/\mu_0$. Thus $q$ is the usual magnetization surface-charge density.
It is an equivalent source for $\mathbf H$, not a magnetic monopole in $\mathbf B$.

## 3. Boundary integral equation

Set $g=\hat{\mathbf n}\cdot\mathbf H_0$. Normal $\mathbf B$ continuity gives

$$
\mu_m\left(g-K'q-\frac12q\right)
=\mu_0\left(g-K'q+\frac12q\right),
$$

and hence

$$
\left[\frac{\mu_m+\mu_0}{2}I+(\mu_m-\mu_0)K'\right]q
=(\mu_m-\mu_0)g.
$$

Define

$$
\mu_r=\frac{\mu_m}{\mu_0},\qquad
\beta=\frac{\mu_r-1}{\mu_r+1}.
$$

The normalized second-kind equation is

$$
\boxed{(I+2\beta K')q=2\beta g.}
$$

Equivalently, when $U_0$ is available, $g=-\partial_n U_0$ and

$$
(I+2\beta K')q=-2\beta\partial_n U_0.
$$

For several disjoint bodies in vacuum, use the $\beta$ belonging to the metal at
each target row. All source surfaces contribute to $K'$. Interfaces shared by
two non-vacuum materials require a generalized interface description and are
deferred.

### Hollow shields

The same equation applies on inner and outer surfaces, with normals always
out of the metal. In block form,

$$
\begin{pmatrix}
I+2\beta K'_{ii} & 2\beta K'_{io}\\
2\beta K'_{oi} & I+2\beta K'_{oo}
\end{pmatrix}
\begin{pmatrix}q_i\\q_o\end{pmatrix}
=2\beta\begin{pmatrix}g_i\\g_o\end{pmatrix}.
$$

Here the first index denotes the target surface and the second the source
surface. The cross-surface blocks retain finite-thickness interactions. An
open-ended shield still bounds a closed metal volume: inner and outer surfaces
must be joined by end/rim surfaces.

### Charge conservation and high permeability

For each connected metal body without magnetic sources,

$$
\int_{\partial\Omega_m}q\,\mathrm dS=0.
$$

This follows from zero net normal flux of $\mathbf H$ on both sides of the complete
metal boundary. Track conservation explicitly in the discretization. A hollow
body's complete boundary includes its cavity surfaces.

At $\mu_r=1$, $\beta=0$ and $q=0$. As $\mu_r$ grows, $\beta$ approaches $1$. For
geometries without free currents linking loops in the metal, bounded-flux
solutions approach $\mathbf H_m=0$ and constant $U$ in each connected metal body,
recovering the ideal-permeability condition.

The limiting operator can have charge modes with vanishing singular values.
Finite but large permeability, and thin shells, can therefore be numerically
demanding despite the second-kind form. Measure conditioning and conservation;
do not assume that bounded $\beta$ ensures uniform conditioning. Treat an exact
infinite-permeability solve separately, with the appropriate charge constraints.
Linked-current topologies require separate treatment of the limiting behavior.

## 4. Initial discretization and field evaluation

Use piecewise constant $q_j$ on each flat triangle $T_j$ and collocate at triangle
centroids $\mathbf x_i$. This gives one unknown per face:

$$
\begin{aligned}
K'_{ij}&=\operatorname{p.v.}\!\int_{T_j}
\partial_{n_i}G(\mathbf x_i,\mathbf y)\,\mathrm dS_{\mathbf y},\\
g_i&=\hat{\mathbf n}_i\cdot\mathbf H_0(\mathbf x_i),\\
A_{ij}&=\delta_{ij}+2\beta_i K'_{ij},\\
\mathbf A\mathbf q&=2\operatorname{diag}(\boldsymbol\beta)\mathbf g.
\end{aligned}
$$

For one common permeability, the right-hand side reduces to $2\beta\mathbf g$.

The integral includes the source triangle area. Do not multiply it by area a
second time. For a planar triangle collocated at its interior centroid,
$K'_{ii}=0$; the jump term is already represented by the identity in $\mathbf A$. These
statements do not prescribe vertex-collocation solid-angle coefficients.

Assemble potential and gradient evaluation operators away from the surface:

$$
\begin{aligned}
S_{pj}&=\int_{T_j}G(\mathbf x_p,\mathbf y)\,\mathrm dS_{\mathbf y},\\
D_{paj}&=\int_{T_j}\partial_{x_a}G(\mathbf x_p,\mathbf y)\,\mathrm dS_{\mathbf y},\\
u_p&=\sum_j S_{pj}q_j,\\
H_{pa}&=H_{0,pa}-\sum_j D_{paj}q_j,\\
B_{pa}&=\mu(\mathbf x_p)H_{pa}.
\end{aligned}
$$

Proposed array conventions: $S$ has shape `(Npoints, Nfaces)`, $D$ has shape
`(Npoints, 3, Nfaces)`, and $q$ has shape `(Nfaces,)` or `(Nfaces, Nrhs)`.
The boundary operator is the target-normal contraction of the principal-value
gradient operator, with its self term handled explicitly.

On-surface evaluation needs an explicit metal/vacuum trace selection and the
corresponding jump. Do not approximate traces using an arbitrary fixed offset.
At edges and corners the smooth-surface half-jump formula does not directly
define a unique pointwise field; initially restrict boundary trace evaluation
to face interiors.

Use analytic triangle integrals where practical. Near interactions between
closely spaced shell surfaces require stable analytic evaluation or controlled
adaptive quadrature. A centroid point-source approximation is not sufficient
for these interactions. Begin with dense assembly and direct solves; retain
target chunking to limit temporary memory.

## 5. Relationship to existing bfieldtools code

- `examples/publication_physics/mumetal_shield_example.py` solves an ideal
  shielding problem using an equivalent stream function and scalar-potential
  collocation. It is already a BEM, based on a double-layer representation.
- `bfieldtools/mesh_magnetics.py::scalar_potential_coupling` maps vertex stream
  functions through dipole/double-layer kernels. It is not the single-layer
  operator required here. In the paper, Eqs. (1) and (6) relate the stream
  function to the surface current and double-layer potential as

  $$
  \mathbf j=\nabla_{\parallel}\psi\times\hat{\mathbf n},\qquad
  U(\mathbf x)=\int_S\psi(\mathbf y)
  \hat{\mathbf n}(\mathbf y)\cdot\nabla_{\mathbf y}G(\mathbf x,\mathbf y)
  \,\mathrm dS_{\mathbf y}.
  $$

  The new single-layer density $q$ is distinct from the stream function $\psi$.
- `bfieldtools/integrals.py::triangle_potential_uniform` evaluates the integral
  of $1/r$ over a triangle and is a starting point for $S$. Apply the $1/(4\pi)$
  normalization explicitly.
- Audit `gamma0`, `omega`, and related existing kernels for an analytic
  gradient implementation. Establish signs and normalizations by independent
  numerical checks before reuse.
- Existing coil $\mathbf B$ coupling supplies $\mathbf B_0$ in tesla. Use
  $\mathbf H_0=\mathbf B_0/\mu_0$ as the excitation. With
  `scalar_potential_coupling(..., multiply_coeff=False)`, the scalar coupling
  has the paper's $1/(4\pi)$ normalization; `multiply_coeff=True` instead
  includes $\mu_0/(4\pi)$. Verify that the existing coupling obeys
  $\mathbf B_0=-\mu_0\nabla U_0$ when reusing it.

## 6. Staged implementation plan

### Stage 1: geometry and integral operators (implemented)

1. Add `bfieldtools/bem/__init__.py` and a focused operator module, provisionally
   `operators.py`. Keep the new face-density basis distinct from vertex stream
   functions.
2. Update package configuration: `pyproject.toml` currently explicitly lists
   only `bfieldtools`; ensure the new `bfieldtools.bem` package is included in
   built distributions.
3. Define a boundary description carrying face ordering, normals, connected
   metal-body labels, and relative permeability per body.
4. Validate watertightness, nondegenerate faces, and consistent orientation.
   Preserve metal-outward normals on cavity surfaces; independently orienting
   every closed surface outward from its enclosed volume would be wrong.
   Detect or reject intersections and unsupported shared-material interfaces.
5. Implement single-layer potential and gradient coupling, including the
   $1/(4\pi)$ normalization and target chunking.
6. Assemble $K'$ at face centroids with the explicit zero self term. Validate
   close-surface interactions before accepting thin-shell geometries.

### Stage 2: finite-permeability solver

1. Add a solver module, provisionally `magnetostatics.py`, that accepts a
   validated boundary, positive finite $\mu_r$ values, and $\mathbf H_0$ at the centroids.
2. Assemble $\mathbf A$ and solve for one or multiple excitation columns. Support reuse
   of the factorization for fixed geometry and permeability.
3. Return $q$ together with residuals and area-weighted charge errors per metal
   body. Assess a conservative or constrained discretization if collocation
   charge drift materially affects large-$\mu_r$ accuracy; do not silently alter
   $q$ after solving.
4. Add potential and $\mathbf H$ evaluation from the induced density. Define $\mathbf B$ evaluation
   with explicit region/permeability information before adding automatic point
   classification, including cavities and boundary-side handling.
5. Document the paper's negative-gradient potential convention and units in the public API.
   Finalize public names after the minimal solver and evaluation path works.

### Stage 3: validation and numerical limits

1. Check analytic triangle kernels against independent quadrature away from
   singularities. Check $\nabla S$ against finite differences at off-surface
   points and verify the normal-derivative jump from both sides.
2. Verify $\mu_r=1$ gives exactly zero response, plus linearity and consistent
   multiple-right-hand-side results.
3. For a solid sphere of radius $a$ in a uniform $\mathbf H_0$, verify the interior field

   $$
   \mathbf H_{\mathrm{inside}}=\frac{3}{\mu_r+2}\mathbf H_0
   $$

   and the exterior dipole moment

   $$
   \mathbf m=4\pi a^3\frac{\mu_r-1}{\mu_r+2}\mathbf H_0.
   $$

   With the present convention, the exact surface density is

   $$
   q=3\frac{\mu_r-1}{\mu_r+2}
   \bigl(\mathbf H_0\cdot\hat{\mathbf n}\bigr).
   $$

   Measure convergence with mesh refinement, including field samples away
   from the boundary and the exterior dipole coefficient.
4. Derive the concentric spherical-shell solution by matching the $l=1$
   harmonic potentials at both radii. Validate cavity attenuation, metal
   fields, and the exterior dipole over permeability and thickness sweeps.
   This exercises both boundary components and their normal orientations.
5. Check tangential $\mathbf H$ and normal $\mathbf B$ continuity using independent one-sided
   evaluations, rather than only the assembled-system residual.
6. Check total charge per metal body and far-field decay; test interactions
   between multiple bodies. Check covariance under rigid transformations.
7. Compare large-$\mu_r$ cavity fields with the existing ideal-shield calculation
   on compatible geometries. Report conditioning and conservation errors as
   functions of mesh size, shell thickness, and $\mu_r$.

Use these physical and integral-identity tests as acceptance criteria. Set
numerical tolerances from measured convergence; solving the discrete linear
system accurately alone does not establish physical accuracy.

### Stage 4: examples and integration

1. Add a spherical-shell validation example and a finite-thickness cylindrical
   shield example showing finite-permeability versus ideal shielding.
2. Connect coil stream functions through their existing $\mathbf B_0$ coupling, converting
   to $\mathbf H_0$ for the solve. Expose the induced response as a reusable linear map
   so it can participate in coil optimization.
3. Document geometry requirements, units, supported material assumptions,
   trace conventions, and measured thin-shell/high-permeability limitations.
4. Run focused BEM and reused-kernel tests, then relevant existing integral and
   magnetic coupling tests. Check that the subpackage is present in a built
   wheel and imports successfully from that wheel.

## 7. Later extensions

Consider Galerkin discretization with explicit mass matrices and conservation,
higher-order densities, matrix-free iterative solves and acceleration, shared
interfaces between general piecewise-constant media, and a separately derived
thin-shell model. Nonlinear permeability requires additional modeling and
iteration; the present uniform-permeability linear equation is not sufficient.

## References

- A. J. Mäkinen et al., *Magnetic-field modeling with surface currents. Part I.
  Physical and computational principles of bfieldtools*, Journal of Applied
  Physics 128, 063906 (2020).
  [Published article](https://doi.org/10.1063/5.0016090);
  [open preprint](https://arxiv.org/html/2005.10060).
  Sec. II B and Eqs. (1), (4), and (6) establish the scalar-potential and
  stream-function conventions used here. The finite-permeability single-layer
  equation in this plan is a derivation extending those conventions.
- E. Y. Tsymbal, *Section 6: Magnetostatics*: macroscopic field equations,
  interface conditions, scalar potential, and magnetization surface sources.
  https://tsymbal.unl.edu/sites/unl.edu.cas.physics.tsymbal/files/media/file/section6-Magnetostatics.pdf
- A. S. Ferguson, X. Zhang, and G. Stroink, *A complete linear discretization
  for calculating the magnetic field using the boundary element method*,
  IEEE Transactions on Biomedical Engineering 41(5), 455–460 (1994).
  https://doi.org/10.1109/10.293220
  Referenced by the existing uniform-triangle potential implementation.

The single-layer signs and normalized equation above follow the paper's
$\mathbf H=-\nabla U$ convention and the metal-outward-normal convention; translate
conventions before comparing formulas from other sources.
