# Volume-current formulation

Physical interpretation, compatible current spaces, and evolution equations.
For implementation status and remaining work, see [PLAN.md](PLAN.md).
Local normalization and matrix assembly are derived in
[EDGE_CURRENT_MATRICES.md](EDGE_CURRENT_MATRICES.md); geometry-only
integration is described in [INDUCTANCE_INTEGRALS.md](INDUCTANCE_INTEGRALS.md).

## Equivalent magnetization: the surface-to-volume connection

The current vector potential is an equivalent magnetic dipole moment density:

$$
\mathbf M_{\mathrm{eq}}=\mathbf T,
\qquad [\mathbf T]=\mathrm{A/m},
\qquad d\mathbf m=\mathbf T\,dV.
$$

This extends a construction already used by bfieldtools. With the surface
convention in `mesh_calculus.gradient_matrix(rotated=True)`,

$$
\mathbf K=\nabla_s\psi\times\mathbf n,
\qquad d\mathbf m=\psi\mathbf n\,dS.
$$

Thus the stream function is a normal-oriented dipole moment per unit area
(units A); the volume potential is a vector dipole moment per unit volume
(units A/m). A thin layer with normal magnetization integrates to the surface
density: $\int\mathbf T\,dn=\psi\mathbf n$. The surface stream function is
constant on each connected boundary component. This enforces zero current
across the surface boundary; different boundary components may have different
constant values.

The connection is used numerically in
[`mesh_magnetics.scalar_potential_coupling`](../mesh_magnetics.py): it maps the
linear stream function to a dipole-layer scalar potential using
`triangle_potential_dipole_linear`, with optional far-field dipole approximations
from [`integrals.py`](../integrals.py). The volume assembly currently uses cell
currents rather than dipole kernels; the interpretation explains their physical
equivalence without claiming a separate magnetization API exists.

For a magnetization supported in $\Omega$, the equivalent currents are

$$
\mathbf J=\nabla\times\mathbf T\quad\text{in }\Omega,
\qquad
\mathbf K_{\partial\Omega}=\mathbf T\times\mathbf n.
$$

These are the standard [volume and surface magnetization currents](https://farside.ph.utexas.edu/teaching/jk1/lectures/node60.html).
The implemented `tree_cotree_gauge` sets every boundary-edge coefficient to
zero. Edge integrals determine the tangential trace of the lowest-order
Nédélec field on each triangular face, so this enforces
$\mathbf T_h\times\mathbf n=0$, not merely zero normal current.
Tangential continuity across internal faces also eliminates internal current
sheets. Extending $\mathbf T_h$ by zero outside the conductor therefore gives
exactly the represented volume current as its distributional curl, with no
additional sheet currents. The normal component of $\mathbf T_h$ may jump;
this does not create a sheet current.

Boundary conditions and gauge fixing are separate steps in the same
construction: boundary coefficients vanish first, then the interior tree
coefficients vanish. The retained cotree coefficients uniquely determine
$\mathbf T_h$ and $\mathbf J_h$ in the supported solid-ball components.
There is no remaining gradient freedom in this discrete representation.
Equivalent magnetization here describes the current; it is not a constitutive
magnetization of the conductor. The physical model still assumes permeability
$\mu_0$.

Two useful consequences of the zero tangential trace are

$$
\mathbf m=\frac12\int_\Omega\mathbf r\times\mathbf J_h\,dV
=\int_\Omega\mathbf T_h\,dV,
$$

$$
\int_\Omega\mathbf J_h\cdot\mathbf A_{\mathrm{ext}}\,dV
=\int_\Omega\mathbf T_h\cdot\mathbf B_{\mathrm{ext}}\,dV,
\qquad \mathbf B_{\mathrm{ext}}=\nabla\times\mathbf A_{\mathrm{ext}}.
$$

The second identity follows by integration by parts: the boundary term
vanishes. In particular, the implemented loop-coupling integral can also be
understood as coupling an equivalent dipole density to the loop's magnetic
field. Pointwise reconstruction of $\mathbf T_h$, a dipole-moment operator,
and evaluation of this alternative coupling integral are not yet exposed.

## Exterior magnetic scalar potential (not yet implemented)

The equivalent magnetization gives a scalar potential outside the conductor:

$$
\boxed{
U(\mathbf r)=\frac{1}{4\pi}\int_\Omega
\mathbf T_h(\mathbf r')\cdot
\frac{\mathbf r-\mathbf r'}{|\mathbf r-\mathbf r'|^3}\,dV',
\qquad \mathbf r\notin\overline\Omega.
}
$$

Choose the additive constant so that $U\to0$ at infinity. With coordinates in
metres and current coefficients in amperes, $U$ has units A and

$$
\mathbf H=-\nabla U,\qquad \mathbf B=-\mu_0\nabla U
\quad\text{outside }\Omega.
$$

This matches the default normalization of the surface
`scalar_potential_coupling`; its `multiply_coeff=True` convention instead
returns $\mu_0 U$. The formula describes the conductor's current contribution;
any externally imposed field must be added separately. The supported
solid-ball components admit a single-valued exterior potential. Extending to
handle-spanning currents requires revisiting global scalar-potential existence
and possible cuts, alongside the topology extensions.

The zero tangential boundary trace and internal tangential continuity established
above ensure that this dipole distribution represents exactly the implemented
volume current. It does not introduce additional current sheets. The exterior
relation $\mathbf B=-\mu_0\nabla U$ must not be used inside current-carrying
cells: the equivalent-magnetization field there also contains a local
$\mu_0\mathbf T_h$ term.

### Exact cell reduction to existing triangle integrals

Each lowest-order Whitney/Nédélec field on a tetrahedron $K$ has the form

$$
\mathbf T_K(\mathbf r')=\mathbf a_K+\mathbf b_K\times\mathbf r',
\qquad \nabla\times\mathbf T_K=2\mathbf b_K.
$$

For example, with affine barycentric functions $\lambda_i$ and constant
gradients $\mathbf g_i=\nabla\lambda_i$, the edge function is
$\mathbf N_{ij}=\lambda_i\mathbf g_j-\lambda_j\mathbf g_i$ and its rotational
part has $\mathbf b_{ij}=\mathbf g_i\times\mathbf g_j$.
Consequently,

$$
\bigl[\mathbf T_K(\mathbf r')-\mathbf T_K(\mathbf r)\bigr]
\cdot(\mathbf r-\mathbf r')=0.
$$

Here $\mathbf T_K(\mathbf r)$ is the affine continuation of the cell field to
the observation point, even though that point is outside $K$. It is not the
physical value of the zero-extended magnetization at that exterior point.

Define the uniform triangle potential and its cell face sum by

$$
S_f(\mathbf r)=\int_f\frac{dS'}{|\mathbf r-\mathbf r'|},
\qquad
\mathbf Q_K(\mathbf r)=\sum_{f\subset\partial K}\mathbf n_{Kf}S_f(\mathbf r),
$$

where $\mathbf n_{Kf}$ is the unit outward normal of face $f$ relative to $K$.
Since $\nabla_{\mathbf r'}|\mathbf r-\mathbf r'|^{-1}
=(\mathbf r-\mathbf r')/|\mathbf r-\mathbf r'|^3$, the divergence theorem gives

$$
\mathbf Q_K(\mathbf r)=\int_K
\frac{\mathbf r-\mathbf r'}{|\mathbf r-\mathbf r'|^3}\,dV'.
$$

The scalar potential therefore reduces to

$$
\boxed{
U(\mathbf r)=\frac{1}{4\pi}\sum_K
\mathbf T_K(\mathbf r)\cdot\mathbf Q_K(\mathbf r).
}
$$

The existing [`triangle_potential_uniform`](../integrals.py) evaluates $S_f$.
Thus this is an analytic source reduction requiring no volume quadrature.
Internal faces must remain in the cell sum: the normal component of
$\mathbf T_h$ may jump across them, so their contributions do not generally
cancel. The boundary condition sets the tangential, not the normal, component
of $\mathbf T_h$ to zero.

### Proposed coupling operator and validation

For global oriented edge functions, one could assemble

$$
C_{pe}=\frac{1}{4\pi}\sum_{K\supset e}
\mathbf N_e|_K(\mathbf r_p)\cdot\mathbf Q_K(\mathbf r_p).
$$

With the existing coefficient prolongation $t=Px$, exterior values would then
be $U=CPx$. Local edge functions must include their local-to-global orientation
signs. The tree–cotree prolongation already zeros boundary and tree edges.
Construction needs geometry and the chosen edge-potential representation;
`CellCurrentBasis` alone does not retain the gauge or potential functions.

This reduction is documented here as a proposed operator, not an available
volume API. Before implementing it, establish:

- Independent volume-dipole quadrature checks, including near-boundary points.
- Agreement of $-\mu_0\nabla U$ with independent Biot–Savart evaluation of the
  reconstructed currents, using exterior finite differences or analytic gradients.
- Far-field agreement with the total dipole moment $\int_\Omega\mathbf T_h\,dV$.
- Orientation, geometric scaling, and invariance of exterior values under
  compatible changes of gauge representing the same physical current.
- Stable coordinate scaling and a controlled far-field treatment. Face sums
  and affine extrapolation can suffer cancellation far from a cell; the
  analytic identity alone is not a floating-point accuracy guarantee. See the
  related stabilization in [INDUCTANCE_INTEGRALS.md](INDUCTANCE_INTEGRALS.md#implemented-numerical-evaluator).

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

The implemented solver unknowns are **gauge-fixed independent edge coefficients**.
For each solid-ball conductor, the implementation uses a boundary-compatible **tree gauge**:
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

The tests validate the implemented gauge-fixed edge basis against the full
conserved face-flux space using dense linear algebra. An explicit dual-cycle
basis and global-pattern construction are not implemented. These are useful
reference constructions for extending beyond solid-ball components.

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

The implemented evaluator integrates the finite source potential over testing
cells and averages source/testing directions. It handles self and touching
pairs without coincident-point kernel evaluation or diagonal adjustments;
see [INDUCTANCE_INTEGRALS.md](INDUCTANCE_INTEGRALS.md) for its convergence
estimator and limitations. General volume-current vector-potential and
magnetic-field evaluation operators remain to be implemented.

## Modes and time evolution

After removing redundant coefficients, solve

$$
Rv_k=\lambda_kMv_k, \qquad v_i^\mathsf TMv_j=\delta_{ij}.
$$

For the implemented passive model, eigenvalues are positive decay rates, with
$\tau_k=1/\lambda_k$. The examples use SciPy's dense symmetric generalized
eigensolver on small problems.

For an externally imposed vector potential,

$$
M\dot x+Rx=f(t), \qquad
f_a(t)=-\int_\Omega\mathbf j_a\cdot\partial_t\mathbf A_{\mathrm{ext}}\,dV.
$$

With $x=\sum_k v_k a_k$, modal evolution is

$$
\dot a_k+\lambda_k a_k=v_k^\mathsf Tf(t).
$$

For the implemented prescribed-loop drive, write
$\mathbf A_{\mathrm{ext}}=I(t)\mathbf A_{\mathrm{loop}}^{(1)}$ and define

$$
m_a=\int_\Omega\mathbf j_a\cdot\mathbf A_{\mathrm{loop}}^{(1)}\,dV,
\qquad
M\dot x+Rx=-m\dot I(t).
$$

`loop_coupling` returns $m$ for a closed oriented polygon outside the conductor.
`modal_response` integrates
$\dot a_k+\lambda_k a_k=-(v_k^{\mathsf T}m)\dot I(t)$, with optional initial
amplitudes at the first requested time. A constant source current supplies no
continuing excitation in this stationary-conductor model. The source waveform
is prescribed; source-circuit back reaction is not included.

Check reduced simulations against direct integration of the full discrete system.
For free decay, verify the energy identity

$$
\frac{d}{dt}\left(\tfrac12x^\mathsf TMx\right)=-x^\mathsf TRx.
$$

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

