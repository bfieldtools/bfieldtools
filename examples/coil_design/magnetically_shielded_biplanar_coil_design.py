'''
Magnetically shielded  coil
===========================
Compact example of design of a biplanar coil within a cylindrical shield.
The effect of the shield is prospectively taken into account while designing the coil.
The coil is positioned close to the end of the shield to demonstrate the effect
'''


import numpy as np
from mayavi import mlab
import trimesh


from bfieldtools.mesh_class import MeshWrapper
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops, plot_data_on_vertices

import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


#Load simple plane mesh that is centered on the origin
planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

planemesh.apply_scale(scaling_factor)

#Specify coil plane geometry
center_offset = np.array([9, 0, 0]) * scaling_factor
standoff = np.array([0, 4, 0]) * scaling_factor

#Create coil plane pairs
coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                         planemesh.faces, process=False)

coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                     planemesh.faces, process=False)

joined_planes = coil_plus.union(coil_minus)

#Create mesh class object
coil = MeshWrapper(mesh_obj=joined_planes, fix_normals=True)

# Separate object for shield geometry
shieldmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/closed_cylinder.stl'), process=True)
shieldmesh.apply_scale(15)

shield = MeshWrapper(mesh_obj=shieldmesh, process=True, fix_normals=True)


###############################################################
# Set up target  points and plot geometry

#Here, the target points are on a volumetric grid within a sphere
# Set up target and stray field points

#Here, the target points are on a volumetric grid within a sphere

center = np.array([9, 0, 0]) * scaling_factor

sidelength = 3 * scaling_factor
n = 12
xx = np.linspace(-sidelength/2, sidelength/2, n)
yy = np.linspace(-sidelength/2, sidelength/2, n)
zz = np.linspace(-sidelength/2, sidelength/2, n)
X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

x = X.ravel()
y = Y.ravel()
z = Z.ravel()

target_points = np.array([x, y, z]).T

#Turn cube into sphere by rejecting points "in the corners"
target_points = target_points[np.linalg.norm(target_points, axis=1) < sidelength/2]  + center


#Plot coil, shield and target points

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                size=(800, 800))

coil.plot_mesh(representation='surface')
shield.plot_mesh()
mlab.points3d(*target_points.T)

################################################################
# Let's design a coil without taking the magnetic shield into account

#The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function
target_field = np.zeros(target_points.shape)
target_field[:, 1] = target_field[:, 1] + 1 # Homogeneous Z-field


target_rel_error = np.zeros_like(target_field)
target_rel_error[:, 0] += 0.01

target_abs_error = np.zeros_like(target_field)
target_abs_error[:, 0] += 0.001
target_abs_error[:, 1:3] += 0.005

target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}

import mosek

coil.j, coil.prob = optimize_streamfunctions(coil,
                                   [target_spec],
                                   objective='minimum_inductive_energy',
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                   )



##############################################################
# Plot coil windings and target points

loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=10)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f)

B_target = coil.B_coupling(target_points) @ coil.j

mlab.quiver3d(*target_points.T, *B_target.T)

#################################################################
# Now, let's compute the effect of the shield on the field produced by the coil

# Points slightly inside the shield
d = np.mean(np.diff(shield.mesh.vertices[shield.mesh.faces[:,0:2]],axis=1), axis=0)/10
points = shield.mesh.vertices - d*shield.mesh.vertex_normals

# Calculate primary potential matrix at the shield surface
P_prim = coil.U_coupling(points)

# Calculate linear collocation BEM matrix
P_bem = shield.U_coupling(points)

# Recalculate diag elements according to de Munck paper
#for diag_index in range(P_bem.shape[0]):
#    P_bem[diag_index, diag_index] = 0
#    P_bem[diag_index, diag_index] = -P_bem[diag_index, :].sum()

# Matrix misses one rank, make it invertible
# by rank-one update (sets potential of constant dipole layer)
#P_bem += np.ones(P_bem.shape)/P_bem.shape[0]


# Solve equivalent stream function for the perfect linear mu-metal layer.
# This is the equivalent surface current in the shield that would cause its
# scalar magnetic potential to be constant
shield.j =  np.linalg.solve(P_bem, P_prim @ coil.j)

##########################################################
# Plot the difference in field when taking the shield into account

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

B_target = coil.B_coupling(target_points) @ coil.j

B_target_w_shield = coil.B_coupling(target_points) @ coil.j + shield.B_coupling(target_points) @ shield.j

B_quiver = mlab.quiver3d(*target_points.T, *(B_target_w_shield - B_target).T, colormap='viridis', mode='arrow')
f.scene.isometric_view()
mlab.colorbar(B_quiver, title='Difference in magnetic field (a.u.)')

import seaborn as sns
import matplotlib.pyplot as plt




fig, axes = plt.subplots(1, 3, figsize=(10, 4))

fig.suptitle('Component-wise effect of magnetic shield on target field amplitude distribution')
for ax_idx, ax in enumerate(axes):

    sns.distplot(B_target[:, ax_idx], label='Without shield', ax=ax)
    sns.distplot(B_target_w_shield[:, ax_idx], label='With shield', ax=ax)
    ax.set_xlabel('Magnetic field (a.u.)')

    if ax_idx == 2:
        ax.legend()

fig.tight_layout(rect=[0, 0.03, 1, 0.95])


###############################################################
# Let's redesign the coil taking the shield into account prospectively

shield.coupling = np.linalg.solve(P_bem, P_prim)

secondary_B_coupling = shield.B_coupling(target_points) @ shield.coupling

total_B_coupling = coil.B_coupling(target_points) + secondary_B_coupling

target_spec_w_shield = {'coupling':total_B_coupling, 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}


coil.j2, coil.prob2 = optimize_streamfunctions(coil,
                                   [target_spec_w_shield],
                                   objective='minimum_inductive_energy',
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                   )

##############################################################
# Plot the newly designed coil windings and field at the target points

loops, loop_values= scalar_contour(coil.mesh, coil.j2, N_contours=10)
f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f)

B_target2 = total_B_coupling @ coil.j2
mlab.quiver3d(*target_points.T, *B_target2.T)

###############################################################
# Plot the difference in stream functions

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_data_on_vertices(coil.mesh, np.nan_to_num(100 * (coil.j-coil.j2)/coil.j), figure=f, colorbar=True)

mlab.colorbar(title='Relative error (%)')
