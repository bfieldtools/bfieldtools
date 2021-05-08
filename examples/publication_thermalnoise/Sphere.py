import numpy as np
import matplotlib.pyplot as plt
import trimesh
from mayavi import mlab
from bfieldtools.thermal_noise import (
    compute_current_modes,
    noise_covar,
    noise_var,
    visualize_current_modes,
)
from bfieldtools.mesh_magnetics import magnetic_field_coupling

import pkg_resources


font = {"family": "normal", "weight": "normal", "size": 16}
plt.rc("font", **font)

# Fix the simulation parameters
d = 1e-3
sigma = 3.8e7
T = 293
kB = 1.38064852e-23
mu0 = 4 * np.pi * 1e-7
freqs = np.array((0,))


Nchunks = 8
quad_degree = 2

##############################################################################
# Unit sphere
# ------------


Np = 10
radius = np.linspace(0.1, 1, Np)
fp = np.zeros((1, 3))

B = np.zeros((Np, 3))
for i in range(Np):
    mesh = trimesh.load(
        pkg_resources.resource_filename("bfieldtools", "example_meshes/unit_sphere.stl")
    )
    mesh.apply_scale(radius[i])

    B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)

    vl = compute_current_modes(
        obj=mesh,
        T=T,
        resistivity=1 / sigma,
        thickness=d,
        mode="AC",
        return_eigenvals=False,
    )

    #    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
    #               size=(800, 800))
    #    visualize_current_modes(mesh,vl[:,:,0], 8, 1)

    vl[:, 0] = np.zeros(vl[:, 0].shape)  # fix DC-component

    Btemp = noise_var(B_coupling, vl)

    B[i] = Btemp[:, :, 0]

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
scene.scene.z_minus_view()
surface = scene.children[0].children[0].children[0].children[0]
surface.actor.property.representation = "wireframe"
surface.actor.mapper.scalar_visibility = False
scene.scene.camera.position = [0.0, 0.0, -5.530686305704514]
scene.scene.camera.focal_point = [0.0, 0.0, 0.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 1.0, 0.0]
scene.scene.camera.clipping_range = [3.485379442647469, 8.118646600290083]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
scene.scene.camera.position = [0.0, 0.0, -4.570815128681416]
scene.scene.camera.focal_point = [0.0, 0.0, 0.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 1.0, 0.0]
scene.scene.camera.clipping_range = [2.535106977394602, 7.1443773556116374]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()


Ban = mu0 * np.sqrt(2 * sigma * d * kB * T / (3 * np.pi * (radius) ** 2))

plt.figure(figsize=(5, 5))
plt.semilogy(radius, Ban * 1e15, linewidth=2, label="Analytic")
plt.semilogy(
    radius,
    np.sqrt(B[:, 2]) * 1e15,
    "x",
    markersize=10,
    markeredgewidth=2,
    label="Numerical",
)
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False)
plt.xlabel("R (m)")
plt.ylabel(r"$B_z$ noise at DC (fT/rHz)")
plt.ylim(5, 100)
plt.tight_layout()


RE = np.abs((np.sqrt(B[:, 2]) - Ban)) / np.abs(Ban) * 100
plt.figure()
plt.plot(radius, RE)
plt.xlabel("Sphere radius")
plt.ylabel("Relative error (%)")
