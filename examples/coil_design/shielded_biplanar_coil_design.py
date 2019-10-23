'''
Shielded biplanar coil design
=============================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes. In addition, the coils have an outer surface
for which (in a linear fashion) a secondary current is created, which zeroes the
normal component of the field produced by the primary coil at the secondary coil
surface. The combination of the primary and secondary coil currents are specified to create
the target field, and their combined inductive energy is minimized.

NB. The secondary coil current is entirely a function of the primary coil current
and the geometry.
'''

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import trimesh

from bfieldtools.mesh_class import MeshWrapper, CouplingMatrix
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


#Load simple plane mesh that is centered on the origin
planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane.obj'), process=False)

planemesh.apply_scale(scaling_factor*1.6)

#Specify coil plane geometry
center_offset = np.array([0, 0, 0]) * scaling_factor
standoff = np.array([0, 5, 0]) * scaling_factor

#Create coil plane pairs
coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                         planemesh.faces, process=False)

coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                     planemesh.faces, process=False)

joined_planes = coil_plus.union(coil_minus)

#Create mesh class object
coil = MeshWrapper(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)

shieldmesh = joined_planes.copy()
shieldmesh.vertices *= np.array([1.5, 1.5, 1.5])

shieldcoil = MeshWrapper(verts=shieldmesh.vertices, tris=shieldmesh.faces, fix_normals=True)



##############################################################
# Compute inductances and coupling


M11 = coil.inductance
M22 = shieldcoil.inductance
# Constrain boundary to zero and consider only inneverts
M11 = M11#[coil.inner_verts][:, coil.inner_verts]
M22 = M22[shieldcoil.inner_verts][:, shieldcoil.inner_verts]
# Add rank-one matrix, so that M22 can be inverted (for zero mean functions)
#M22 += np.ones_like(M22)/M22.shape[0]
#M11 += np.ones_like(M11)/M11.shape[0]

from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix_from_A

M21 = mutual_inductance_matrix_from_A(shieldcoil.mesh, coil.mesh)
M21 = M21[shieldcoil.inner_verts]

# Mapping from I1 to I2, constraining flux through shieldcoil to zero
P = -np.linalg.solve(M22, M21)



##############################################################
# Set up target and stray field points

#Here, the target points are on a volumetric grid within a sphere

center = np.array([0, 0, 0]) * scaling_factor

sidelength = 2 * scaling_factor
n = 8
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


##############################################################
# Compute C matrices that are used to compute the generated magnetic field

coil.C = CouplingMatrix(coil, compute_C)
shieldcoil.C = CouplingMatrix(shieldcoil, compute_C)

##############################################################
# Create bfield specifications used when optimizing the coil geometry

#The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function

target_field = np.zeros(target_points.shape)
target_field[:, 1] = target_field[:, 1] + 1

target_rel_error = np.zeros_like(target_field)
target_rel_error[:, 0] += 0.01

target_abs_error = np.zeros_like(target_field)
target_abs_error[:, 0] += 0.001
target_abs_error[:, 1:3] += 0.005

target_spec = {'C':coil.C(target_points) + shieldcoil.C(target_points)[:, :, shieldcoil.inner_verts]@P, 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target_field':target_field}
#[:, :, coil.inner_verts]

objective_matrix = M11 - M21.T @ np.linalg.pinv(M22) @ M21

##############################################################
# Run QP solver
import mosek

coil.I, prob = optimize_streamfunctions(coil,
                                   [target_spec],
                                   objective=objective_matrix,
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}},
                                   boundary_constraints='all_zero'
                                   )

shieldcoil.I = np.zeros((len(shieldcoil.mesh.vertices, )))

shieldcoil.I[shieldcoil.inner_verts] = P @ coil.I

from bfieldtools.viz import plot_data_on_vertices

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))

plot_data_on_vertices(coil.mesh, coil.I, figure=f)
plot_data_on_vertices(shieldcoil.mesh, shieldcoil.I, figure=f)

#############################################################
# Plot coil windings and target points

N_contours = 10

loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=N_contours)
sloops, sloop_values= scalar_contour(shieldcoil.mesh, shieldcoil.I, N_contours=N_contours)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f)
plot_3d_current_loops(sloops, colors='auto', figure=f)

B_target = coil.C(target_points) @ coil.I + shieldcoil.C(target_points) @shieldcoil.I

mlab.quiver3d(*target_points.T, *B_target.T)




extent = 30




x1 = np.linspace(-extent, extent, 101) * scaling_factor

y1 = z1 = np.zeros_like(x1)

line1_points = np.vstack((x1, y1, z1)).T

B_line1 = coil.C(line1_points) @ coil.I + shieldcoil.C(line1_points) @ shieldcoil.I


y2 = np.linspace(-extent, extent, 101) * scaling_factor

z2 = x2 = np.zeros_like(y2)

line2_points = np.vstack((x2, y2, z2)).T

B_line2 = coil.C(line2_points) @ coil.I + shieldcoil.C(line2_points) @ shieldcoil.I



z3 = np.linspace(-extent, extent, 101) * scaling_factor

x3 = y3 = np.zeros_like(z1)

line3_points = np.vstack((x3, y3, z3)).T


B_line3 = coil.C(line3_points) @ coil.I + shieldcoil.C(line3_points) @ shieldcoil.I

fig, axes = plt.subplots(1, 1)

for ax_idx, ax in enumerate([axes]):
    ax.semilogy(x1 / scaling_factor, np.linalg.norm(B_line1, axis=-1), label='X')
    ax.semilogy(y2 / scaling_factor, np.linalg.norm(B_line2, axis=-1), label='Y')
    ax.semilogy(z3 / scaling_factor, np.linalg.norm(B_line3, axis=-1), label='Z')
    ax.set_title('Field component %d'% ax_idx)

plt.ylabel('Field amplitude (target field units)')
plt.xlabel('Distance from origin')
plt.grid(True, which='minor', axis='y')
plt.grid(True, which='major', axis='y', color='k')
plt.grid(True, which='major', axis='x')

plt.legend()


plt.show()

#%%
from bfieldtools.magnetic_field_mesh import compute_C_analytic, compute_U

x = y = np.linspace(-20, 20, 50)
X,Y = np.meshgrid(x, y, indexing='ij')
points = np.zeros((X.flatten().shape[0], 3))
points[:, 0] = X.flatten()
points[:, 1] = Y.flatten()

CB1 = compute_C_analytic(coil.mesh, points)
CB2 = compute_C_analytic(shieldcoil.mesh, points)

CU1 = compute_U(coil.mesh, points)
CU2 = compute_U(shieldcoil.mesh, points)

B1 = CB1 @ coil.I
B2 = CB2 @ shieldcoil.I

U1 = CU1 @ coil.I
U2 = CU2 @ shieldcoil.I


#%% Plot
B = (B1.T + B2.T)[:2].reshape(2, x.shape[0], y.shape[0])
lw = np.sqrt(B[0]**2 + B[1]**2)
lw = 2*lw/np.max(lw)
xx = np.linspace(-1,1, 16)
#seed_points = 0.51*np.array([xx, -np.sqrt(1-xx**2)])
#seed_points = np.hstack([seed_points, (0.51*np.array([xx, np.sqrt(1-xx**2)]))])
#plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
#               start_points=seed_points.T, integration_direction='both')
U = (U1 + U2).reshape(x.shape[0], y.shape[0])
U /= np.max(U)
plt.figure()
plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
           extent=(x.min(), x.max(), y.min(), y.max()))
plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
               #start_points=seed_points.T,
               integration_direction='both')

cc1 = scalar_contour(coil.mesh, coil.mesh.vertices[:,2], contours= [-0.001])[0][0]
cc2 = scalar_contour(shieldcoil.mesh, shieldcoil.mesh.vertices[:,2], contours= [-0.001])[0][0]

plt.plot(cc1[:,1], cc1[:,0], linewidth=3.0)
plt.plot(cc2[:,1], cc2[:,0], linewidth=3.0)

plt.xticks([])
plt.yticks([])



