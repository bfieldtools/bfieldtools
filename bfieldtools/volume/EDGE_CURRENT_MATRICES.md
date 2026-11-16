# Resistance and inductance from tetrahedral edge curls

Status: mathematical draft for discussion; no solver implementation yet.
This develops the current-vector-potential formulation proposed in
[PLAN.md](PLAN.md): lowest-order, first-kind Nédélec edge functions, their
piecewise-constant curls, and a boundary-compatible tree–cotree gauge.
Here “first order” means the six-function Whitney edge element, denoted
`Nedelec_0` in the plan. All lengths below use SI units when computing physical
resistance and inductance.

## Local basis and normalization

Let a tetrahedron have vertices $\mathbf r_0,\ldots,\mathbf r_3$, ordered so that

$$
6V_T=\det[\mathbf r_1-\mathbf r_0,\;
\mathbf r_2-\mathbf r_0,\;\mathbf r_3-\mathbf r_0]>0.
$$

Write $\lambda_i$ for its barycentric coordinates and
$\mathbf g_i=\nabla\lambda_i$. For the directed edge $i\to j$, use

$$
\mathbf N_{ij}=\lambda_i\mathbf g_j-\lambda_j\mathbf g_i,
\qquad
\int_{i\to j}\mathbf N_{ij}\cdot d\mathbf l=1.
$$

The integral along every other local edge is zero. Consequently the
coefficient $t_{ij}$ is an edge integral of the current vector potential and
has units A. There is **no edge-length multiplier** in this normalization.
The basis current is

$$
\boxed{\mathbf c_{T,ij}=\nabla\times\mathbf N_{ij}
=2\mathbf g_i\times\mathbf g_j.}
$$

It is constant within $T$ and has units $\mathrm{m}^{-2}$, so
$t_{ij}\mathbf c_{T,ij}$ has units $\mathrm{A/m^2}$.

## Expressing the curl with the opposing edge

Choose exactly the directed local edge list

$$
(e_0,e_1,e_2,e_3,e_4,e_5)=(01,02,03,23,31,12),
\qquad
\mathbf d_{T,n}=\mathbf r_j-\mathbf r_i\quad(e_n=ij).
$$

The opposing-edge index is

$$
p(n)=(n+3)\bmod6.
$$

The parentheses matter in code: `(n + 3) % 6`.
For four distinct indices, the geometric identity is

$$
\mathbf g_i\times\mathbf g_j
=\varepsilon_{ijkl}\frac{\mathbf r_l-\mathbf r_k}{6V_T},
\qquad \varepsilon_{0123}=+1.
$$

Reversing $13$ to $31$ in the edge list makes all three complementary pairs
$(01,23)$, $(02,31)$, and $(03,12)$ positively oriented. Thus the directed
opposing edge gives the current direction without a local sign factor:

$$
\boxed{
\mathbf c_{T,n}=\frac{\mathbf d_{T,p(n)}}{3V_T}.
}
$$

| $n$ | Potential edge | Opposing directed edge | $3V_T\,\mathbf c_{T,n}$ |
| --- | --- | --- | --- |
| 0 | $0\to1$ | $2\to3$ | $\mathbf r_3-\mathbf r_2$ |
| 1 | $0\to2$ | $3\to1$ | $\mathbf r_1-\mathbf r_3$ |
| 2 | $0\to3$ | $1\to2$ | $\mathbf r_2-\mathbf r_1$ |
| 3 | $2\to3$ | $0\to1$ | $\mathbf r_1-\mathbf r_0$ |
| 4 | $3\to1$ | $0\to2$ | $\mathbf r_2-\mathbf r_0$ |
| 5 | $1\to2$ | $0\to3$ | $\mathbf r_3-\mathbf r_0$ |

One way to see the direction is that $\lambda_i$ and $\lambda_j$ are both
constant along the opposing edge. Their gradients are perpendicular to it,
so their cross product is parallel to it. The determinant fixes the sign and
the factor $1/(6V_T)$.

For example, on $(\mathbf r_0,\mathbf r_1,\mathbf r_2,\mathbf r_3)
=(\mathbf0,\mathbf e_x,\mathbf e_y,\mathbf e_z)$, $V_T=1/6$ and
$\mathbf c_{T,01}=2(0,-1,1)$, whereas
$\mathbf c_{T,02}=2(1,0,-1)$, along the directed opposing edge $3\to1$.

## Global orientation and assembly

Assign each mesh edge one global direction, for example from the smaller global
vertex index to the larger. Let $e(T,n)$ be the global edge corresponding to
local edge $n$, and $s_{T,n}=+1$ if their directions agree, $-1$ otherwise.
The global basis current associated with edge $e$ is

$$
\mathbf j_e|_T=
\begin{cases}
s_{T,n}\mathbf c_{T,n},&e=e(T,n),\\
\mathbf0,&e\not\subset T.
\end{cases}
$$

The sign $s_{T,n}$ belongs to the **potential edge**, not the opposing edge.
The opposing-edge vector in the formula remains in the local direction in the
table. Replacing that vector with a globally directed one would require its
own orientation conversion.

For compact formulas define

$$
\mathbf b_{T,e}=
\begin{cases}
s_{T,n}\mathbf d_{T,p(n)},&e=e(T,n),\\
\mathbf0,&e\not\subset T,
\end{cases}
\qquad
\mathbf j_e|_T=\frac{\mathbf b_{T,e}}{3V_T}.
$$

The current reconstruction is
$\mathbf J_h|_T=\sum_e t_e\mathbf b_{T,e}/(3V_T)$.
Normal components match across faces because these currents are curls of
globally conforming edge functions. Each cell current is divergence-free.

## Resistance matrix

Let $\rho_T=1/\sigma_T>0$ be scalar resistivity, constant in each cell. The
Joule loss is $P=t^{\mathsf T}R^{\mathrm{edge}}t$, where

$$
R^{\mathrm{edge}}_{ef}
=\int_\Omega\rho\,\mathbf j_e\cdot\mathbf j_f\,dV
=\sum_T\rho_TV_T\,\mathbf j_e|_T\cdot\mathbf j_f|_T.
$$

Substituting the opposing-edge vectors gives the exact entry

$$
\boxed{
R^{\mathrm{edge}}_{ef}
=\frac19\sum_T\frac{\rho_T}{V_T}
\mathbf b_{T,e}\cdot\mathbf b_{T,f}.
}
$$

Equivalently, the local $6\times6$ matrix before global orientation signs is

$$
\boxed{
R^T_{nm}=\frac{\rho_T}{9V_T}
\mathbf d_{T,p(n)}\cdot\mathbf d_{T,p(m)}.
}
$$

Scatter it to the global entries with factor $s_{T,n}s_{T,m}$.
Only edges sharing a tetrahedron couple, so resistance is sparse.
For a constant symmetric resistivity tensor, replace
$\rho_T\,\mathbf b_{T,e}\cdot\mathbf b_{T,f}$ by
$\mathbf b_{T,e}^{\mathsf T}\boldsymbol\rho_T\mathbf b_{T,f}$.

## Inductance matrix

For nonmagnetic conductors in free space, define the scalar cell-pair integral

$$
K_{TU}=\int_T\int_U\frac{dV_{\mathbf y}\,dV_{\mathbf x}}
{|\mathbf x-\mathbf y|}.
$$

The magnetic energy is $W=\tfrac12t^{\mathsf T}M^{\mathrm{edge}}t$, with

$$
M^{\mathrm{edge}}_{ef}=\frac{\mu_0}{4\pi}
\sum_{T,U}K_{TU}\,\mathbf j_e|_T\cdot\mathbf j_f|_U.
$$

Hence

$$
\boxed{
M^{\mathrm{edge}}_{ef}
=\frac{\mu_0}{36\pi}
\sum_{T,U}\frac{K_{TU}}{V_TV_U}
\mathbf b_{T,e}\cdot\mathbf b_{U,f}.
}
$$

The corresponding local $6\times6$ interaction block is

$$
\boxed{
M^{TU}_{nm}=\frac{\mu_0K_{TU}}{36\pi V_TV_U}
\mathbf d_{T,p(n)}\cdot\mathbf d_{U,p(m)}.
}
$$

Scatter with $s_{T,n}s_{U,m}$. The sum is over all ordered cell pairs,
including $T=U$; there is no extra factor $1/2$ in a matrix entry.
The $1/2$ belongs to the energy. If evaluating only unordered distinct pairs,
assemble both transposed blocks.

The scalar $K_{TU}$ is exactly the geometric integral discussed in
[INDUCTANCE_INTEGRALS.md](INDUCTANCE_INTEGRALS.md). That note gives analytic
source integration and an exact reduction to face-pair distance integrals;
it does not yet establish a fully integrated formula for arbitrary cell pairs.
The vector dot products above are exact. Numerical integration is confined to
the scalar cell interactions. Inductance is generally dense.

## Boundary conditions and the tree–cotree gauge

For this first construction, assume a connected conductor homeomorphic to a
solid ball, with a conforming mesh. Set all boundary-edge coefficients to zero.
This sets the tangential trace of $\mathbf T_h$ to zero and therefore enforces
zero normal current on boundary faces. The remaining candidate unknowns are
edges not contained in the boundary surface, including interior edges whose
endpoints both happen to be boundary vertices.

These coefficients still contain gradient freedom:

$$
t_{ij}\longmapsto t_{ij}+\psi_j-\psi_i.
$$

To preserve the zero boundary coefficients, take $\psi$ constant on the
connected boundary and choose that constant as zero. A concrete compatible
tree–cotree construction is:

1. Contract all boundary vertices to a single root in the mesh edge graph.
2. Choose a spanning tree of this contracted graph, using the remaining edges.
   Edges that become self-loops cannot belong to the tree.
3. Set the tree-edge coefficients to zero by choosing the interior nodal
   values of $\psi$. Retain all other nonboundary edges as cotree unknowns.

There is one tree edge per interior vertex. For this topology the retained
edge curls form an independent basis for the conserved, zero-boundary-flux
current space. Thus a retained basis function is simply its original global
edge curl; the gauge does not require adding tree-edge currents to it.
The tree supplies the return path for interpreting its coefficient as a
circulation of $\mathbf T_h$.

More generally, write $t=Px$ for any chosen independent representation after
boundary constraints and gauge fixing. Then

$$
\boxed{R=P^{\mathsf T}R^{\mathrm{edge}}P,\qquad
M=P^{\mathsf T}M^{\mathrm{edge}}P.}
$$

For the simple tree gauge, $P$ just selects the retained edges, so one can
assemble their submatrices directly. With cavities or handles, the gradient
kernel, boundary conditions, and possible missing global current patterns
must be treated explicitly as described in [PLAN.md](PLAN.md); the solid-ball
construction above is not a general topology algorithm.

## Computational precedents and reference-cell equivalence

[Tong, Lombeyda, Hirani and Desbrun (2003), *Discrete Multiscale Vector Field
Decomposition*, Appendix C](https://geometry.caltech.edu/pubs/TLHD03.pdf)
explicitly derives the opposing-edge identity
$V_T\nabla\lambda_i\times\nabla\lambda_j=\mathbf d_{kl}/6$, with the
opposing edge oriented along the cross product. They use it in their discrete
calculus on tetrahedral meshes. Their vector-potential discretization differs
from the Whitney edge potential here, so this is a precedent for the geometric
identity, not for this entire eddy-current solver.

[Passarotto, Pitassi and Specogna (2022), *Foundations of volume integral methods
for eddy current problems*, Section 2.4.1, equations (40)–(45)](https://www.comphys.com/pdf/cmame_22_vinco.pdf)
describes the CARIDDI computational formulation: edge-integral potential
coefficients, currents obtained as curls of Nédélec functions, resistance and
inductance integrals, and tree–cotree gauging. This supports the computational
formulation; the cited section does not specify our opposing-edge implementation
or our local edge ordering.

The direct geometric formula is also equivalent to the usual affine
reference-cell construction. To see this, let
$\mathbf x=\mathbf r_0+B\widehat{\mathbf x}$, with $\det B=6V_T>0$ for the
unit reference tetrahedron. The covariant transformation and its curl give

$$
\mathbf N(\mathbf x)=B^{-\mathsf T}\widehat{\mathbf N}(\widehat{\mathbf x}),
\qquad
\nabla\times\mathbf N
=\frac{B}{\det B}\bigl(\widehat\nabla\times\widehat{\mathbf N}\bigr).
$$

With our edge list, the reference curl is $2\widehat{\mathbf d}_{p(n)}$.
Since $B\widehat{\mathbf d}_{p(n)}=\mathbf d_{T,p(n)}$, the transformed curl
is exactly $\mathbf d_{T,p(n)}/(3V_T)$. This derivation explains why the
reference-cell machinery can be collapsed into edge vectors and cell volumes
for these affine, lowest-order tetrahedra. Curved cells do not have the same
constant-curl simplification.

## Checks before implementation

- Compare every local curl $2\mathbf g_i\times\mathbf g_j$ with the directed
  opposing-edge formula, including permuted positive cell orderings.
- Verify oriented face fluxes of each edge curl against the face–edge incidence
  matrix; their magnitudes on faces containing the edge are one.
- Check zero boundary flux and independence after the tree gauge, with retained
  dimension $N_{\mathrm{interior\ edges}}-N_{\mathrm{interior\ vertices}}$.
- Check symmetry and positive definiteness of the reduced matrices. The full
  edge matrices contain a curl kernel and are only positive semidefinite.
- Check $R$ has units $\Omega$ and $M$ has units H. Under uniform geometric
  scaling by $a>0$ at fixed resistivity and fixed coefficient normalization,
  $R\mapsto R/a$ and $M\mapsto aM$.

The key convention to settle is therefore unit edge-integral normalization
with the directed local edge list above. It gives the constant current
$\mathbf d_{T,p(n)}/(3V_T)$ and fixes all factors in both matrices.
