"""
Field of a Helmholtz coils pair
==========================================================

Example on how to compute the magnetic field from current line segments forming a Helmholtz coil pair.

Visualization of the 3D field using Mayavi
"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import polyline, structured_field

pv.set_plot_theme("document")


import matplotlib.pyplot as plt
import numpy as np

from bfieldtools.utils import cylinder_points
from bfieldtools.line_magnetics import magnetic_field

# Create helmholtz coil with radius R
R = 5

c_points = cylinder_points(
    radius=R, length=0, nlength=1, nalpha=100, orientation=np.array([0, 1, 0])
)
c_points[:, 1] = 0
c_points = np.vstack((c_points, c_points[0, :]))


c1_points = c_points - np.array([0, R / 2, 0])
c2_points = c_points + np.array([0, R / 2, 0])

plotter = pv.Plotter()
plotter.add_mesh(polyline(*c1_points.T).tube(radius=0.01), show_scalar_bar=False)
plotter.add_mesh(polyline(*c2_points.T).tube(radius=0.01), show_scalar_bar=False)

box = 3 * R
n = 50

xx = np.linspace(-box, box, n)
yy = np.linspace(-box, box, n)
zz = np.linspace(-box, box, n)
X, Y, Z = np.meshgrid(xx, yy, zz, indexing="ij")

x = X.ravel()
y = Y.ravel()
z = Z.ravel()

b_points = np.array([x, y, z]).T

B = np.zeros(b_points.shape)

B += magnetic_field(c1_points, b_points)
B += magnetic_field(c2_points, b_points)


B_matrix = B.reshape((n, n, n, 3))

B_matrix_norm = np.linalg.norm(B_matrix, axis=-1)

#%% Visualization using PyVista

plotter = pv.Plotter()

field = structured_field(
    X,
    Y,
    Z,
    B_matrix[:, :, :, 0],
    B_matrix[:, :, :, 1],
    B_matrix[:, :, :, 2],
    scalars=B_matrix_norm,
)

vectors = plotter.add_mesh(
    field.glyph(orient="vectors", scale="magnitude", factor=X[1, 0, 0] - X[0, 0, 0]),
    show_scalar_bar=False,
)


vcp = plotter.add_mesh(
    field.slice(normal="y").glyph(orient="vectors", scale="magnitude"),
    show_scalar_bar=False,
)
# For prettier picture:

iso = plotter.add_mesh(
    field.contour(10, scalars="magnitude"),
    opacity=0.2,
    cmap="viridis",
    show_scalar_bar=False,
)


# A trick to make transparency look better: cull the front face
mesh_actor(iso).prop.frontface_culling = True

# Settings


#%% Plot field magnitude on x and y axes using matplotlib

plt.figure()

z1 = np.linspace(0, 30, 31)

x1 = y1 = np.zeros_like(z1)

line1_points = np.vstack((x1, y1, z1)).T


Bh_line1 = magnetic_field(c1_points, line1_points) + magnetic_field(
    c2_points, line1_points
)

plt.semilogy(
    z1,
    np.linalg.norm(Bh_line1, axis=1) / np.linalg.norm(Bh_line1, axis=1)[0],
    label="Helmholtz, z-axis",
)


y2 = np.linspace(0, 30, 31)

z2 = x2 = np.zeros_like(y2)

line2_points = np.vstack((x2, y2, z2)).T

Bh_line2 = magnetic_field(c1_points, line2_points) + magnetic_field(
    c2_points, line2_points
)


plt.semilogy(
    y2,
    np.linalg.norm(Bh_line2, axis=1) / np.linalg.norm(Bh_line2, axis=1)[0],
    label="Helmholtz, y-axis",
)
plt.ylabel("Field amplitude (target field units)")
plt.xlabel("Distance from origin")
plt.grid(True, which="minor", axis="y")
plt.grid(True, which="major", axis="y", color="k")
plt.grid(True, which="major", axis="x")

plt.legend()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
