"""

Compare SUH and SPH basis functions for the magnetic field
=============================================================

"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import triangle_mesh, vector_glyphs

pv.set_plot_theme("document")


import numpy as np
import matplotlib.pyplot as plt
import trimesh
from bfieldtools import resources


from bfieldtools.mesh_magnetics import magnetic_field_coupling
from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic
from bfieldtools.mesh_magnetics import scalar_potential_coupling
from bfieldtools.sphtools import compute_sphcoeffs_mesh, basis_fields
from bfieldtools.suhtools import SuhBasis
from bfieldtools.utils import load_example_mesh


mesh = load_example_mesh("bunny_repaired")

mesh.vertices -= mesh.vertices.mean(axis=0)

mesh_field = mesh.copy()
mesh_field.vertices += 0.005 * mesh_field.vertex_normals
mesh_field = trimesh.smoothing.filter_laplacian(mesh_field, iterations=1)

Ca, Cb = basis_fields(mesh_field.vertices, 4)

bsuh = SuhBasis(mesh, 25)
Csuh = magnetic_field_coupling_analytic(mesh, mesh_field.vertices) @ bsuh.basis

#%%
def plot_basis_fields(C, comps):
    d = 0.17
    i = 0
    j = 0
    for n in comps:
        p = 1.05 * mesh_field.vertices.copy()
        p2 = mesh_field.vertices.copy()
        p[:, 0] += i * d
        p2[:, 0] += i * d
        m = np.max(np.linalg.norm(C[:, :, n], axis=0))
        vectors = plotter.add_mesh(
            vector_glyphs(*p.T, *C[:, :, n].T, mode="arrow"),
            cmap="Greys",
            clim=[0, m],
            show_scalar_bar=False,
        )
        m = np.max(abs((C[:, :, n].T * mesh_field.vertex_normals.T).sum(axis=0)))
        s = plotter.add_mesh(
            triangle_mesh(*p2.T, mesh.faces),
            scalars=(C[:, :, n].T * mesh_field.vertex_normals.T).sum(axis=0),
            cmap="bwr",
            clim=[-m, m],
            show_scalar_bar=False,
        )
        mesh_actor(s).mapper.interpolate_before_map = True
        mesh_actor(s).mapper.lookup_table.n_values = 15
        i += 1


comps = [0, 4, 10, 15]
plotter = scene = pv.Plotter(window_size=(1200, 350))
plot_basis_fields(Ca, comps)
scene.parallel_projection = True
scene.view_xy()
scene.camera.zoom(4)


plotter = scene = pv.Plotter(window_size=(1200, 350))
plot_basis_fields(Csuh, comps)
scene.parallel_projection = True
scene.view_xy()
scene.camera.zoom(4)


#%% test representation
# from bfieldtools.mesh_magnetics import scalar_potential_coupling
# scaling_factor = 0.02
##Load simple plane mesh that is centered on the origin
# file_obj=resources.files('bfieldtools', 'example_meshes/10x10_plane_hires.obj')
# plane = trimesh.load(file_obj=file_obj, process=False)
# plane.apply_scale(scaling_factor)
## Rotate to x-plane
# t = np.eye(4)
# t[1:3,1:3] = np.array([[0,1],[-1,0]])
# plane.apply_transform(t)
# plane = plane.subdivide()
#
# u = scalar_potential_coupling(mesh, plane.vertices, multiply_coeff=False)
# mask = 1 + np.sum(u, axis=1)/(4*np.pi)
# mask[mask < 1e-6]  = 0
#
#
# Ca, Cb = basis_fields(plane.vertices, 6)
# bsuh = SuhBasis(mesh, 48)
# CB = magnetic_field_coupling_analytic(mesh, plane.vertices)
# Csuh = CB @ bsuh.basis
#
# A,B = compute_sphcoeffs_mesh(mesh, 6)
# b1 = np.einsum('ij,lik->lkj', A , Ca) # Mapping from verts to alpha to field
# b2 = np.einsum('ij,klj->lki', bsuh.mass @ bsuh.basis  , Csuh) # Mapping from verts to suh to field
#
# ind=1932
# scalars = np.zeros(mesh.vertices.shape[0])
# scalars[ind] = 1
#                     scalars = mask*np.linalg.norm(CB[:,:,ind], axis=1))
#
#                     scalars = mask*np.linalg.norm(b1[:,:,ind], axis=0))
#
#                     scalars = mask*np.linalg.norm(b2[:,:,ind], axis=0))

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
