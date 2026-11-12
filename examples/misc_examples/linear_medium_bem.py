"""
Solving linearly polarizing material in a homogenenous field
------------------------------------------------------------

The example is based on the following integral equation of
a piecewise harmonic potential :math:`U` on the boundary of the polarizing material

.. math::
    \\frac{\\mu_0 + \\mu}{2} U(\\vec{r}) = \\mu_0 U_\\infty(\\vec{r}) + \\frac{\\mu_0 - \\mu}{4\\pi}\\int_{S} U(\\vec{r}\\,') (\\hat{n}' \\cdot \\nabla') \\frac{1}{|\\vec{r}-\\vec{r}\\,'|} dS'\\,.

* :math:`U` is the magnetic scalar potential to be solved
* :math:`\\mu_0` is the permeability of the free space
* :math:`\\mu` is the permeability of the material
* :math:`U_\\infty` is the magnetic scalar potential without the material
* :math:`S` is the surface of the material

After discretization :math:`U(\\vec{r})=\sum_n u_n h_n(\\vec{r})` and integrating the integral equation over every hat :math:`h_n(\\vec{r})`,
a linear system of equations is obtained

.. math::
    \\mathbf{N}\\mathbf{u} = \\frac{2}{1+\\mu_\\mathrm{r}} \\mathbf{u}_\\infty + \\frac{1 - \\mu_\\mathrm{r}}{4\\pi(1 + \\mu_\\mathrm{r})}\\mathbf{D}\\mathbf{u}\,.

* :math:`\\mathbf{u}` is :math:`U` at the mesh nodes
* :math:`\\mathbf{u}_\\infty` is the potential :math:`U_\\infty` integrated over hat functions
* :math:`\\mu_\\mathrm{r}` is the relative permeability of the material
* :math:`\\mathbf{D}` is the double-layer coupling matrix
* :math:`\\mathbf{N}` is the mesh mass matrix
"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import triangle_mesh, vector_glyphs, structured_field

pv.set_plot_theme("document")


import numpy as np
import matplotlib.pyplot as plt

from bfieldtools.mesh_conductor import MeshConductor as Conductor
from bfieldtools.mesh_magnetics import magnetic_field_coupling
from bfieldtools.mesh_magnetics import scalar_potential_coupling
from bfieldtools.mesh_calculus import mass_matrix
from bfieldtools.utils import load_example_mesh


from trimesh.creation import icosphere

# Use a sphere
# mesh = icosphere(3, 1)

# Use a cube
mesh = load_example_mesh("cube")
mesh.vertices -= mesh.vertices.mean(axis=0)
mesh.apply_scale(0.2)


def Dmatrix(mesh1, mesh2, Nchunks=100):
    """
    "Double-layer potential" coupling between two meshes
    using a Galerkin method with hat basis

    Discretize integral equations using hat functions
    on on both meshes. Potential from mesh1 hat functions
    is calculated analytically and integrated over the
    hat functions of mesh2 numerically

    Parameters
    ----------
    mesh1 : Trimesh object
    mesh2 : Trimesh object
    Nchunks : int, optional
        Number of chunks in the potential calculation. The default is 100.

    Returns
    -------
    None.

    """
    face_points = mesh2.vertices[mesh2.faces]
    weights = np.array([[0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5]])
    # Combine vertices for quadrature points
    Rquad = np.einsum("...ij,ik->...kj", face_points, weights)
    R = Rquad.reshape(-1, 3)
    U = scalar_potential_coupling(
        mesh1, R, Nchunks, multiply_coeff=False, approx_far=True
    )

    face_areas = mesh2.area_faces
    # Reshape and multiply by quadrature weights
    Dcomps = U.reshape(Rquad.shape[:2] + (len(mesh2.vertices),)) * (
        face_areas[:, None, None] / 3
    )
    # Sum up the quadrature points
    D = mesh2.faces_sparse @ Dcomps[:, 0, :]
    D += mesh2.faces_sparse @ Dcomps[:, 1, :]
    D += mesh2.faces_sparse @ Dcomps[:, 2, :]

    if mesh1 is mesh2:
        # Recalculate diagonals
        d = np.diag(D)
        D -= np.diag(d)
        # Make rows sum to -2*pi*(vertex area), should be more accurate
        d2 = -2 * np.pi * mass_matrix(mesh2, lumped=False) - np.diag(D.sum(axis=1))
        D += d2
        # Make D solvable by adding rank-1 matrix
        D += np.ones_like(D) * np.max(np.linalg.svd(D, False, False)) / D.shape[1]

    return D


# Some linear input potentials -> uniform field
def phi0x(r):
    return r[:, 0]


def phi0y(r):
    return r[:, 1]


def phi0z(r):
    return r[:, 2]


def project_to_hats(mesh, func):
    """
    Numerically integrate func over hat functions

    Parameters
    ----------
    mesh : Trimesh object
        the domain for hat functions
    func : function
        potential function for phi0, takes (N,3) array
        of points as input

    Returns
    -------
    p_mat : (Nvertices,) array
        func projected on the hat functions

    """
    # Index vertex points for each faces
    face_points = mesh.vertices[mesh.faces]
    weights = np.array([[0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5]])
    # Combine vertices for quadrature points
    Rquad = np.einsum("...ij,ik->...kj", face_points, weights)
    R = Rquad.reshape(-1, 3)

    # Evaluation func at quadrature points
    p = func(R)
    face_areas = mesh.area_faces
    # Reshape and multiply by quadrature weights
    pcomps = p.reshape(Rquad.shape[:2]) * (face_areas[:, None] / 3)
    # Sum up the quadrature points
    p_mat = mesh.faces_sparse @ pcomps[:, 0]
    p_mat += mesh.faces_sparse @ pcomps[:, 1]
    p_mat += mesh.faces_sparse @ pcomps[:, 2]

    return p_mat


#%% Compute matrices
print("Computing D matrix")
D = Dmatrix(mesh, mesh)
print("Computing mass matrix")
M = mass_matrix(mesh)
Ml = mass_matrix(mesh, lumped=True)


#%% Solve linear BEM equations
print("Computing input potential")
pp = project_to_hats(mesh, phi0x)
mu_r = 100
c1 = (1 - mu_r) / (1 + mu_r)
u = np.linalg.solve(M - 1 / (2 * np.pi) * c1 * D, -2 * pp / (mu_r + 1))

# Plot potential on the mesh
plotter = pv.Plotter()
m = plotter.add_mesh(
    triangle_mesh(*mesh.vertices.T, mesh.faces),
    scalars=u,
    cmap="bwr",
    show_scalar_bar=False,
)
mesh_actor(m).mapper.interpolate_before_map = True
mesh_actor(m).mapper.lookup_table.n_values = 32
#%% Compute fields on a plane
# Load plane for the visualization of the potential
plane = load_example_mesh("10x10_plane_hires", process=True)
t = np.eye(4)
t[1:3, 1:3] = np.array([[0, 1], [-1, 0]])
plane.apply_transform(t)
plane.apply_scale(0.3)
plane = plane.subdivide()

Uplane = scalar_potential_coupling(
    mesh, plane.vertices, multiply_coeff=False, approx_far=True
)

uprim = phi0x(plane.vertices)
usec = (mu_r - 1) / (4 * np.pi) * Uplane @ u
uplane = uprim + usec

# Meshgrid on the same plane for the bfield
X, Y = np.meshgrid(
    np.linspace(-1.5, 1.5, 50), np.linspace(-1.5, 1.5, 50), indexing="ij"
)
pp = np.zeros((50 * 50, 3))
pp[:, 0] = X.flatten()
pp[:, 1] = Y.flatten()

Bplane = magnetic_field_coupling(mesh, pp, analytic=True)
bprim = pp * 0  # copy pp
# add x directional field
mu0 = 1e-7 * 4 * np.pi
bprim[:, 0] = -1
# In this simulation mu0==1 is assumed
# magnetic_field_coupling uses mu0 in SI units
bplane = (mu_r - 1) / (4 * np.pi) / mu0 * Bplane @ u + bprim
# uplane[Uplane.sum(axis=1)]

#%% Plot the results
plotter = pv.Plotter(window_size=(1000, 1000))
m = plotter.add_mesh(
    triangle_mesh(*plane.vertices.T, plane.faces),
    scalars=uplane,
    cmap="bwr",
    show_scalar_bar=False,
)
mesh_actor(m).mapper.interpolate_before_map = True
mesh_actor(m).mapper.lookup_table.n_values = 32
#                         *(-gradient(uplane, plane)), color=(0,0,0),
#               scale_mode='none', scale_factor=0.01, mode='arrow')
plotter.view_xy()

plotter = pv.Plotter(window_size=(1000, 1000))
vectors2 = plotter.add_mesh(
    vector_glyphs(
        *pp.T, *bplane.T, scale_factor=0.02, mode="arrow", scale_mode="vector"
    ),
    color=(1, 0, 0),
    show_scalar_bar=False,
)
plotter.view_xy()

# Streamline plot for the cube
bplane[:, 2] = 0
vecfield = structured_field(
    *pp.T.reshape(3, 50, 50, 1), *bplane.T.reshape(3, 50, 50, 1)
)
vecnorm = vecfield
streams = []

Nq = 40
q = np.zeros((Nq, 3))
q[:, 1] = np.linspace(-1.5, 1.5, Nq)
q[:, 0] = 1.5
extent = np.array([-1.5, 1.5, -1.5, 1.5, 0, 0])
for qi in q:
    stream = plotter.add_mesh(
        vecnorm.streamlines_from_source(
            pv.PolyData(qi.reshape(1, 3)),
            vectors="vectors",
            integration_direction="both",
            terminal_speed=1e-13,
            max_length=6.0,
        ),
        cmap="viridis",
        show_scalar_bar=False,
    )
    streams.append(stream)

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
