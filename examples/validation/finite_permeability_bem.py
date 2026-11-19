"""
Finite-permeability sphere and hollow shield
==========================================

Use the stage-one BEM operators to solve the constant-face collocation equation
in a uniform applied field. Compare a solid sphere and a spherical shell with
their analytic solutions. The potential convention is H = -grad(U).

The dense solve is written explicitly here; this example does not require a
solver class. Refining the flat triangles improves agreement with the smooth
spheres. In particular, coarse meshes are inaccurate at very high permeability.

Run from the repository root with::

    python examples/validation/finite_permeability_bem.py

For a headless run, set MPLBACKEND=Agg. No mesh files or 3D renderer are needed.
"""

import matplotlib.pyplot as plt
import numpy as np
import trimesh

from bfieldtools.bem import (
    MagneticBoundary,
    adjoint_double_layer,
    single_layer_gradient,
    single_layer_potential,
)


def solve_uniform_field(boundary):
    """Solve (I + 2 beta K')q = 2 beta n.H_0 for H_0 = z_hat A/m."""
    k = adjoint_double_layer(boundary)
    g = boundary.mesh.face_normals[:, 2]
    matrix = np.eye(len(g)) + 2 * boundary.beta[:, None] * k
    rhs = 2 * boundary.beta * g
    q = np.linalg.solve(matrix, rhs)
    residual = np.linalg.norm(np.einsum("ij,j->i", matrix, q) - rhs) / np.linalg.norm(
        rhs
    )
    charge = np.sum(boundary.mesh.area_faces * q)
    relative_charge = abs(charge) / np.sum(boundary.mesh.area_faces * abs(q))
    print(
        "%d faces: relative residual %.2e, relative net charge %.2e"
        % (len(q), residual, relative_charge)
    )
    return q


def shell_coefficients(inner_radius, outer_radius, mu_r):
    """Coefficients of U = C*z, A*z+B*z/r^3, -z+c*z/r^3.

    Match U and mu*d_r(U) at both spherical interfaces. The incident H field
    has unit magnitude in the z direction; potentials are in amperes.
    """
    a, b = inner_radius, outer_radius
    matrix = np.array(
        [
            [-a, a, a**-2, 0],
            [-1, mu_r, -2 * mu_r / a**3, 0],
            [0, b, b**-2, -(b**-2)],
            [0, mu_r, -2 * mu_r / b**3, 2 / b**3],
        ]
    )
    return np.linalg.solve(matrix, [0, 0, -b, -1])


# %% Geometry: one connected metal body, with two boundary components for a shell.
mu_r = 10.0
inner_radius = 0.6
outer_radius = 1.0
subdivisions = 3
outer = trimesh.creation.icosphere(subdivisions=subdivisions, radius=outer_radius)
inner = trimesh.creation.icosphere(subdivisions=subdivisions, radius=inner_radius)
inner.invert()  # The inner surface normal must point OUT OF METAL, into the cavity.

solid = MagneticBoundary(outer, mu_r=mu_r)
shell = MagneticBoundary(trimesh.util.concatenate([outer, inner]), mu_r=mu_r)

# Avoid comparing directly at the ideal spherical interfaces: the BEM mesh is
# a faceted approximation, and the two field traces differ there.
z = np.linspace(0, 2.0, 151)
z = z[(abs(z - inner_radius) > 0.035) & (abs(z - outer_radius) > 0.035)]
points = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])

# %% Analytic reference fields and potentials on the positive z axis.
dipole = outer_radius**3 * (mu_r - 1) / (mu_r + 2)
solid_h = np.full_like(z, 3 / (mu_r + 2))
solid_u = -solid_h * z
outside = z > outer_radius
solid_h[outside] = 1 + 2 * dipole / z[outside] ** 3
solid_u[outside] = -z[outside] + dipole / z[outside] ** 2

cavity, metal_a, metal_b, exterior = shell_coefficients(
    inner_radius, outer_radius, mu_r
)
shell_h = np.full_like(z, -cavity)
shell_u = cavity * z
metal = (z > inner_radius) & (z < outer_radius)
shell_h[metal] = -metal_a + 2 * metal_b / z[metal] ** 3
shell_u[metal] = metal_a * z[metal] + metal_b / z[metal] ** 2
shell_h[outside] = 1 + 2 * exterior / z[outside] ** 3
shell_u[outside] = -z[outside] + exterior / z[outside] ** 2

# %% Evaluate H = H_0 - grad(Sq) and U = -z + Sq using the analytic operators.
fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for column, (name, boundary, exact_h, exact_u) in enumerate(
    [
        ("Solid sphere", solid, solid_h, solid_u),
        ("Spherical shell", shell, shell_h, shell_u),
    ]
):
    q = solve_uniform_field(boundary)
    h = np.array([0.0, 0.0, 1.0]) - single_layer_gradient(boundary.mesh, points) @ q
    u = -z + np.einsum("ij,j->i", single_layer_potential(boundary.mesh, points), q)
    # B is mu_0*H in vacuum and mu_0*mu_r*H in the metal. Plot H here so
    # the potential-gradient sign and the interface jump are visible directly.
    print("%s: center H_z = %.6f A/m; analytic %.6f A/m" % (name, h[0, 2], exact_h[0]))
    axes[0, column].plot(z, exact_h, "k--", label="Analytic sphere")
    axes[0, column].plot(z, h[:, 2], label="BEM (flat triangles)")
    axes[1, column].plot(z, exact_u, "k--")
    axes[1, column].plot(z, u)
    axes[0, column].set_title(name + r", $\mu_r=10$")
    axes[0, column].legend()
    for axis in axes[:, column]:
        axis.axvline(outer_radius, color="0.7", linewidth=1)
        if column == 1:
            axis.axvline(inner_radius, color="0.7", linewidth=1)
        axis.grid(alpha=0.3)
    axes[1, column].set_xlabel("z (m)")

axes[0, 0].set_ylabel(r"$H_z$ (A/m)")
axes[1, 0].set_ylabel(r"$U$ (A)")
fig.tight_layout()
plt.show()
