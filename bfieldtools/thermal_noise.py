'''
Contains functions for computing thermal noise in conductive thin objects.

'''

import numpy as np
from scipy.linalg import eigh
from mayavi import mlab

from .mesh_magnetics import magnetic_field_coupling
from .mesh_calculus import laplacian_matrix, mass_matrix, laplacian_matrix_w_holes, mass_matrix_w_holes
from .mesh_inductance import self_inductance_matrix
from . import utils


def compute_current_modes(mesh, boundaries=None):
    '''
    Computes eddy-current modes for a mesh using surface laplacian.
    Uses Dirichlet boundary condition, i.e., stream function is zero at boundary:
    no current flow outside the surface.
    The modes are normalized so that the squared norm of the stream function gradient
    integrates to 1 over the surface. With this normalization, the resistances
    of the current modes are R_k = 1/(sigma*d), sigma = conductivity, d = thickness.
    See Zevenhoven et al. (2014).

    Parameters
    ----------
    mesh: Trimesh mesh object
        The surface mesh
    boundaries: list of N_holes

    Returns
    -------
    vl: Nvertices x Nvertices array
        The normalized eddy-current modes vl[:,i]

    '''
    boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)

    if boundaries:
        L_holes = laplacian_matrix_w_holes(mesh, inner_verts, boundaries)
        M_holes = mass_matrix_w_holes(mesh, inner_verts, boundaries)

        u, v = eigh(-L_holes, M_holes)


        #Normalize the laplacien eigenvectors

        for i in range(v.shape[1]):
            v[:, i] = v[:, i]/np.sqrt(u[i])

        #Assign values per vertex
        vl = np.zeros((mesh.vertices.shape[0], v.shape[1]))

        vl[inner_verts] = v[:-len(boundaries)]

        for b_idx, b in enumerate(boundaries):
            vl[b] = v[len(inner_verts) + b_idx]

    else:
        L = laplacian_matrix(mesh)
        M = mass_matrix(mesh)

        u, v = eigh(-L.todense()[inner_verts][:, inner_verts], M.todense()[inner_verts][:, inner_verts])

        #Normalize the laplacien eigenvectors
        vl = np.zeros(M.shape)
        for i in range(v.shape[1]):
            vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])

    return vl

def compute_dc_Bnoise(mesh, vl, fp, sigma, d, T):
    '''
    Computes the magnetic noise at DC due to thermal motion of charge carriers (Jonhson-Nyquist noise)
    in a relatively thin conductor.

    Parameters
    ----------
    mesh: Trimesh mesh object - the surface mesh
    vl: Nvertices x Nvertices array
        The normalized eddy-current modes vl[:,i]
    fp: Nfieldpoints x 3 array
        Coordinates of the fieldpoints
    sigma: float
        Conductivity of the surface
    d: float
        Thickness of the surface
    T: float
        Temperature of the surface

    Returns
    -------
    B: Nfieldpoints x 3components array
        magnetic RMS noise at DC

    '''

    kB = 1.38064852e-23 #Boltzmann constant

    B_coupling = magnetic_field_coupling(mesh, fp)

    B = np.zeros(fp.shape)
    for i in range(vl.shape[1]):
        vec = vl[:, i] * np.sqrt(4 * kB * T * sigma * d)

        B += (B_coupling @ vec)**2

    B = np.sqrt(B) #RMS

    return B

def compute_dc_Bnoise_covar(mesh, vl, fp, sigma, d, T):
    '''
    Computes the magnetic noise covariance at DC due to thermal motion of charge carriers (Jonhson-Nyquist noise)
    in a relatively thin conductor.

    Parameters
    ----------
    mesh: Trimesh mesh object - the surface mesh
    vl: Nvertices x Nvertices array
        The normalized eddy-current modes vl[:,i]
    fp: Nfieldpoints x 3 array
        Coordinates of the fieldpoints
    sigma: float
        Conductivity of the surface
    d: float
        Thickness of the surface
    T: float
        Temperature of the surface

    Returns
    -------
    B: Nfieldpoints x Nfieldpoints x 3components array
        magnetic noise covariance at DC

    '''

    kB = 1.38064852e-23 #Boltzmann constant

    B_coupling = magnetic_field_coupling(mesh, fp)

    eps = 4*kB*T*sigma*d*np.eye(vl.shape[1])

    B = np.zeros((fp.shape[0], fp.shape[0], 3))
    for i in range(3):
        B[:, :, i] = B_coupling[:, i, :] @ vl @ eps @ (vl.T) @ (B_coupling[:, i, :].T)

    return B


def integrate_Bnoise_covar(B_covar, weighting=None):
    '''
    Computes the (quadrature) integrated noise over a volume spanned by the points in fp.

    Parameters
    ----------
    B_covar: (N_p, N_p, 3) array
        One vector component of the covariance matrix computed by compute_dc_Bnoise_covar
    weighting: (N_p,) array
        Weighting factors for each point in the volume. If None (default), use equal weighting.

    Returns
    -------
    Bnoise_integrated: float
        Integrated noise amplitude over the volume

    '''

    if weighting is None:
        weighting = np.ones((len(B_covar,)))/len(B_covar)
        print('No weighting provided, assuming equal weighting')

    Bnoise_integrated = np.sum(np.outer(weighting, weighting) * B_covar)

    return Bnoise_integrated**0.5


def compute_ac_Bnoise(mesh, vl, fp, freqs, sigma, d, T):
    '''
    Computes the AC magnetic noise due to thermal motion of charge carriers (Jonhson-Nyquist noise)
    in a relatively thin conductor.

    Parameters
    ----------
    mesh: Trimesh mesh object - the surface mesh
    vl: Nvertices x Nvertices array
        The normalized eddy-current modes vl[:,i]
    fp: Nfieldpoints x 3 array
        Coordinates of the fieldpoints
    sigma: float
        Conductivity of the surface
    d: float
        Thickness of the surface
    T: float
        Temperature of the surface

    Returns
    -------
    B: Nfrequencies x Nfieldpoints x 3components array
        magnetic RMS noise across frequencies

    '''

    kB = 1.38064852e-23 #Boltzmann constant

    #Compute field
    B_coupling = magnetic_field_coupling(mesh, fp)

    #Compute mutual inductance at "hat basis"
    Mind = self_inductance_matrix(mesh)
    Mind = 0.5*(Mind+Mind.T) #I think that this enforces M to be symmetric

    #Transform mutual inductance to eddy-current basis
    Mind_lap = vl.T@Mind@vl

    #Eigendecomposition of M
    um_t, vm_t = eigh(Mind_lap)
    um_t = np.flip(um_t)
    vm_t = np.flip(vm_t, axis=1)

    #Let's take just the "99% variance" components to avoid numerical issues
    um = np.zeros(um_t.shape)
    vm = np.zeros(vm_t.shape)

    csum = np.cumsum(um_t**2) / np.sum(um_t**2)
    Ncomps = np.max(np.where(csum < 0.99))

    um[0:Ncomps] = um_t[0:Ncomps]
    vm[0:Ncomps, 0:Ncomps] = vm_t[0:Ncomps, 0:Ncomps]

    #Compute B as a function of frequency
    B = np.zeros((freqs.shape[0], fp.shape[0], fp.shape[1]))

    Rk = 1 / (sigma * d)
    eps = np.sqrt(4*kB*T*Rk)*np.ones((vl.shape[0]))

    for j in range(freqs.shape[0]):
        f = freqs[j]
        omega = 2*np.pi*f


        Rt = Rk/(Rk**2 + omega**2*um**2)
        Rt = np.diag(Rt)

        Lt = um/(Rk**2 + omega**2*um**2)
        Lt = np.diag(Lt)


        currents = -1*(vm@Rt@vm.T@eps + 1j*omega*vm@Lt@vm.T@eps)
        currents = np.abs(currents)

        for i in range(vl.shape[0]):
            vec = currents[i] * vl[:, i]
            B[j, :, 0] += (B_coupling[:, 0, :] @vec)**2
            B[j, :, 1] += (B_coupling[:, 1, :] @ vec)**2
            B[j, :, 2] += (B_coupling[:, 2, :] @ vec)**2
        print("Frequency %f computed" % (f))

    B = np.sqrt(B) #RMS

    return B

def visualize_current_modes(mesh, vl, Nmodes, scale, contours=True, colormap='bwr', dist=0.5):
    '''
    Visualizes current modes up to Nmodes.

    Parameters
    ----------
    mesh: Trimesh mesh object
        The surface mesh
    vl: Nvertices x Nvertices array
        The normalized eddy-current modes vl[:,i]
    Nmodes: int
        Number of modes to be plotted
    scale: float
        Scaling factor
    contours: boolean
        If True, show contours
    colormap: string
        Which (matplotlib) colormap to use

    '''

    N1 = np.floor(np.sqrt(Nmodes))
    dx = (mesh.vertices[:,0].max() - mesh.vertices[:,0].min())*(1+dist)
    dy = (mesh.vertices[:,1].max() - mesh.vertices[:,1].min())*(1+dist)

    i = 0
    j = 0
    for n in range(Nmodes):
        print(i,j)
        points = mesh.vertices.copy()
        points[:,0] += i*dx
        points[:,1] += j*dy
        s = mlab.triangular_mesh(*points.T, mesh.faces,
                             scalars=vl[:, n], colormap=colormap)

        limit = np.max(np.abs(vl[:, n]))

        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.enable_contours = contours

        if i<N1:
            i+=1
        else:
            j+=1
            i=0

    return s
