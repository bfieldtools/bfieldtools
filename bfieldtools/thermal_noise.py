'''
Contains functions for computing thermal noise in conductive thin objects.

'''

import numpy as np
from scipy.linalg import eigh
from mayavi import mlab

from .magnetic_field_mesh import compute_C
from .laplacian_mesh import laplacian_matrix, mass_matrix
from .mutual_inductance_mesh import self_inductance_matrix
from . import utils


def compute_current_modes(mesh):
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

    Returns
    -------
    vl: Nvertices x Nvertices array
        The normalized eddy-current modes vl[:,i]

    '''

    boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)

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

    C = compute_C(mesh, fp)

    B = np.zeros(fp.shape)
    for i in range(vl.shape[1]):
        vec = vl[:, i] * np.sqrt(4 * kB * T * sigma * d)

        B += (C @ vec)**2

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

    C = compute_C(mesh, fp)

    eps = 4*kB*T*sigma*d*np.eye(vl.shape[0])

    B = np.zeros((fp.shape[0], fp.shape[0], 3))
    for i in range(3):
        B[:, :, i] = C[:, :, i] @ vl @ eps @ (vl.T) @ (C[:, :, i].T)

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
    C = compute_C(mesh, fp)

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
            B[j, :, 0] += (C[:, :, 0] @vec)**2
            B[j, :, 1] += (C[:, :, 1] @ vec)**2
            B[j, :, 2] += (C[:, :, 2] @ vec)**2
        print("Frequency %f computed" % (f))

    B = np.sqrt(B) #RMS

    return B

def visualize_current_modes(mesh, vl, Nmodes, scale, contours=True, colormap='bwr'):
    '''
    Visualizes current modes up to Nmodes.
    TODO: make this more flexible.

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
    verts = mesh.vertices
    tris = mesh.faces

    for ii in range(Nmodes):
        n = int(np.sqrt(Nmodes))
        i = ii % n
        j = int(ii/n)

        x = scale*verts[:, 0] + i*(np.max(verts[:, 0]) - np.min(verts[:, 0]))*1.2
        y = scale*verts[:, 1]+ j*(np.max(verts[:, 1]) - np.min(verts[:, 1]))*1.2
        z = scale*verts[:, 2]

        limit = np.max(np.abs(vl[:, ii]))

        s = mlab.triangular_mesh(x, y, z, tris, scalars=vl[:, ii], colormap=colormap)

        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.enable_contours = contours


    return s
