'''
Coil with minimal eddy currents
===============================
Compact example of design of a cylindrical coil surrounded by a RF shield, i.e. a conductive surface.
The effects of eddy currents due to inductive interaction with the shield is minimized
'''
PLOT = False
SAVE_FIGURES = False



import numpy as np
from mayavi import mlab
import trimesh


from bfieldtools.mesh_class import MeshWrapper

from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.mesh_properties import mutual_inductance_matrix_from_A
from bfieldtools.contour import scalar_contour, simplify_contour
from bfieldtools.viz import plot_3d_current_loops, plot_data_on_vertices

import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


#Load example coil mesh that is centered on the origin
coilmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/cylinder.stl'), process=True)

coilmesh.vertices[:,2] -= 3
coilmesh.vertices /= np.max(np.abs(coilmesh.vertices), axis=0)*np.array([2, 2, 1])

angle = np.pi/2
rotation_matrix = np.array([[np.cos(angle), 0, np.sin(angle), 0],
                              [0, 1, 0, 0],
                              [-np.sin(angle), 0, np.cos(angle), 0],
                              [0, 0, 0, 1]
                              ])

coilmesh.apply_transform(rotation_matrix)

coilmesh1 = coilmesh.copy()
#coilmesh1.apply_scale(1.3)

coilmesh2 = coilmesh.copy()

#coilmesh1 = coilmesh.union(coilmesh1)
coilmesh1 = coilmesh1.subdivide()
#coilmesh2 = coilmesh.subdivide()


#Create mesh class object
coil = MeshWrapper(verts=coilmesh1.vertices*0.75, tris=coilmesh1.faces, fix_normals=True)

# Separate object for shield geometry
shield = MeshWrapper(verts=coilmesh2.vertices.copy()*1.1, tris=coilmesh2.faces.copy(), fix_normals=True)
#shield.mesh.vertices[:,2] -= 3
#shield.mesh.vertices *= np.array([1.2, 1.2, 1.2])
#
#angle = np.pi/2
#rotation_matrix = np.array([[np.cos(angle), 0, np.sin(angle), 0],
#                              [0, 1, 0, 0],
#                              [-np.sin(angle), 0, np.cos(angle), 0],
#                              [0, 0, 0, 1]
#                              ])
#
#shield.mesh.apply_transform(rotation_matrix)
#
#shield.mesh = shield.mesh.subdivide()


###############################################################
# Set up target  points and plot geometry

#Here, the target points are on a volumetric grid within a sphere

center = np.array([0, 0, 0])

sidelength = 0.15 * scaling_factor
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
if PLOT:
    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                size=(800, 800))
    coil.plot_mesh()
    shield.plot_mesh()
    mlab.points3d(*target_points.T)




###############################################################
# Compute C matrices that are used to compute the generated magnetic field

mutual_inductance = mutual_inductance_matrix_from_A(coil.mesh, shield.mesh)

# Take into account the field produced by currents induced into the shield
# NB! This expression is for instantaneous step-function switching of coil current, see Eq. 18 in G.N. Peeren, 2003.

shield.M_coupling = np.linalg.solve(-shield.inductance, mutual_inductance.T)
secondary_C = shield.B_coupling(target_points) @ -shield.M_coupling

###############################################################
# Create bfield specifications used when optimizing the coil geometry

#The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function

target_field = np.zeros(target_points.shape)
target_field[:, 1] = target_field[:, 1] + 1

target_abs_error = np.zeros_like(target_field)
target_abs_error[:, 1] += 0.001
target_abs_error[:, 0::2] += 0.005

target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':0, 'abs_error':target_abs_error, 'target':target_field}




def alu_sigma(T):
    ref_T = 293 #K
    ref_rho = 2.82e-8 #ohm*meter
    alpha = 0.0039 #1/K


    rho = alpha * (T - ref_T) * ref_rho + ref_rho

    return 1/rho

resistivity = 1/alu_sigma(T=293) #room-temp Aluminium
thickness = 0.5e-3 # 0.5 mm thick


shield.R = shield.resistance(resistivity=resistivity, thickness=thickness)

#Set boundary vertices to zero
shield.R[shield.boundary_verts, :][:, shield.boundary_verts] = 0

from scipy.linalg import eigh
li, Ui = eigh(shield.R[shield.inner_verts, :][:, shield.inner_verts], shield.inductance[shield.inner_verts, :][:, shield.inner_verts], eigvals=(0, 500))

U = np.zeros((shield.inductance.shape[0], len(li)))
U[shield.inner_verts, :] = Ui


L = np.diag(li)
#
#plt.figure()
#plt.plot(1/li)


#shield.M_coupling = np.linalg.solve(-shield.inductance, mutual_inductance.T)
#secondary_C = shield.B_coupling(target_points) @ -shield.M_coupling


#
#tmin, tmax = 0.001, 0.001
#Fs=10000

time = [0.001, 0.003, 0.005]
eddy_error = [0.01, 0.005, 0.0005]
#time_decay = U @ np.exp(-l[None, :]*time[:, None]) @ np.pinv(U)

time_decay = np.zeros((len(time), shield.inductance.shape[0], shield.inductance.shape[1]))

induction_spec = []


Uinv = np.linalg.pinv(U)
for idx, t in enumerate(time):
     time_decay = U @ np.diag(np.exp(-li*t)) @ Uinv
     eddy_coupling = shield.B_coupling(target_points) @ time_decay @ shield.M_coupling
     induction_spec.append({'coupling':eddy_coupling, 'abs_error':eddy_error[idx], 'rel_error':0, 'target':np.zeros_like(target_field)})

###############################################################
# Run QP solver

import mosek

from bfieldtools.mutual_inductance_mesh import self_inductance_matrix

coil.inductance = self_inductance_matrix(coil.mesh, Nchunks=16, approx=False)

coil.I, prob = optimize_streamfunctions(coil,
                                   [target_spec] + induction_spec,
                                   objective='minimum_inductive_energy',
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                   )

shield.induced_I = shield.M_coupling @ coil.I


###############################################################
# Plot coil windings and target points


loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=20)

#loops = [simplify_contour(loop, min_edge=1e-2, angle_threshold=2e-2, smooth=True) for loop in loops]

#loops = [loop for loop in loops if loop is not None]
if PLOT:
    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.005)

    B_target = coil.B_coupling(target_points) @ coil.I

    mlab.quiver3d(*target_points.T, *B_target.T)

    plot_data_on_vertices(shield.mesh, shield.induced_I, ncolors=256, figure=f, opacity=0.5, cull_back=True)
    plot_data_on_vertices(shield.mesh, shield.induced_I, ncolors=256, figure=f, opacity=1, cull_front=True)

    if SAVE_FIGURES:
        mlab.savefig('/l/bfieldtools/examples/publication_software/eddy_yes.png', figure=f, magnification=4)
#mlab.triangular_mesh(*shield.mesh.vertices.T, shield.mesh.faces, scalars=shield.induced_I)

#mlab.title('Coils which minimize the transient effects of conductive shield')


###############################################################
# For comparison, let's see how the coils look when we ignore the conducting shield


coil.unshielded_I, coil.unshielded_prob = optimize_streamfunctions(coil,
                                   [target_spec],
                                   objective='minimum_inductive_energy',
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                   )

shield.unshielded_induced_I = shield.M_coupling @ coil.unshielded_I

loops, loop_values= scalar_contour(coil.mesh, coil.unshielded_I, N_contours=8)
if PLOT:
    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.005)

    B_target_unshielded = coil.B_coupling(target_points) @ coil.unshielded_I

    mlab.quiver3d(*target_points.T, *B_target_unshielded.T)

    plot_data_on_vertices(shield.mesh, shield.unshielded_induced_I, ncolors=256, figure=f, opacity=0.5, cull_back=True)
    plot_data_on_vertices(shield.mesh, shield.unshielded_induced_I, ncolors=256, figure=f, opacity=1, cull_front=True)

    if SAVE_FIGURES:
        mlab.savefig('/l/bfieldtools/examples/publication_software/eddy_no.png', figure=f, magnification=4)

####################################################################
#Finally, let's compare the time-courses



tmin, tmax = 0, 0.05
Fs=1000

time = np.linspace(tmin, tmax, Fs*(tmax-tmin)+1)

#time_decay = U @ np.exp(-l[None, :]*time[:, None]) @ np.pinv(U)

time_decay = np.zeros((len(time), shield.inductance.shape[0], shield.inductance.shape[1]))

Uinv = np.linalg.pinv(U)
for idx, t in enumerate(time):
     time_decay[idx] = U @ np.diag(np.exp(-li*t)) @ Uinv



B_t = shield.B_coupling(target_points) @ (time_decay @ shield.induced_I).T

unshieldedB_t = shield.B_coupling(target_points) @ (time_decay @ shield.unshielded_induced_I).T

import matplotlib.pyplot as plt

if PLOT and SAVE_FIGURES:
    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(4, 5))

    #ax[0].set_title('Eddy currents minimized')

    axnames = ['X', 'Y', 'Z']
    for i in range(3):

        ax[i].plot(time*1e3, np.mean(B_t, axis=0)[i].T, 'k-', label='Minimized', linewidth=2.0)
        ax[i].plot(time*1e3, np.mean(unshieldedB_t, axis=0)[i].T, 'k--', label='Ignored', linewidth=2.0 )
        #ax[1].set_title('Eddy currents ignored')
    #    ax[i].yaxis.set_major_formatter(ScalarFormatter(useOffset=True))
        #ax[1].set_ylabel('Transient field amplitude')

        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        ax[i].set_ylabel(axnames[i])

        ax[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)

    #ax[2].set_ylabel('Transient field')
    ax[2].set_xlabel('Time (ms)')

    plt.legend()
    fig.tight_layout()
    plt.savefig('/l/bfieldtools/examples/publication_software/eddy_transient_allaxes.pdf')

    plt.close('all')


    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(3, 4))
    ax.plot(time*1e3, np.mean(np.linalg.norm(B_t, axis=1), axis=0).T, 'k-', label='Minimized', linewidth=1.5)
    #ax[0].set_title('Eddy currents minimized')
    ax.set_ylabel('Transient field amplitude')
    ax.plot(time*1e3, np.mean(np.linalg.norm(unshieldedB_t, axis=1), axis=0).T, 'k--', label='Ignored', linewidth=1.5 )
    #ax[1].set_title('Eddy currents ignored')
    ax.set_xlabel('Time (ms)')
    #ax[1].set_ylabel('Transient field amplitude')

    ax.set_ylim(0, 0.15)
    ax.set_xlim(0, 20)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.legend()
    fig.tight_layout()
    plt.savefig('/l/bfieldtools/examples/publication_software/eddy_transient.pdf')


    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(3, 4))
    ax.plot(time*1e3, np.mean(B_t[0], axis=0).T, 'k-', label='Minimized', linewidth=1.5)
    #ax[0].set_title('Eddy currents minimized')
    ax.set_ylabel('Transient field amplitude')
    ax.plot(time*1e3, np.mean(unshieldedB_t[0], axis=0).T, 'k--', label='Ignored', linewidth=1.5 )
    #ax[1].set_title('Eddy currents ignored')
    ax.set_xlabel('Time (ms)')
    #ax[1].set_ylabel('Transient field amplitude')

    ax.set_ylim(-00.05, 0.01)
    ax.set_xlim(0, 20)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.legend()
    fig.tight_layout()
    plt.savefig('/l/bfieldtools/examples/publication_software/eddy_transient_onaxis.pdf')
