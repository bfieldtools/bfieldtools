"""
Functions for constructing real spherical harmonics (Ylms), their gradients
and related magnetic field 'basis vectorfunctions'
(Wlms for r**l components, Vlms for r**(-l) components).

Uses notations and definitions by Plattner and Simons (2014; https://arxiv.org/pdf/1306.3201.pdf)
and the same normalization conventions.

Integration over a surface of unit sphere is
used as the inner product <C,D> = int C dot D dOmega.

Also has many of functions for spherical <-> cartesian transformations.
"""

__all__ = [
    "Blm",
    "Plm",
    "Rotmatrix",
    "SphBasis",
    "Vlm",
    "Wlm",
    "Xlm",
    "basis_fields",
    "basis_potentials",
    "cartesian2spherical",
    "cartvec2sph",
    "compute_sphcoeffs_mesh",
    "derlpmn_em",
    "derxlm",
    "dphiylm",
    "dthylm",
    "field",
    "lpmn_em",
    "plotBVlm_volume",
    "plotBWlm_volume",
    "plotVlm",
    "plotWlm",
    "plotXlm",
    "plotYlm",
    "plotYlms",
    "potential",
    "sinxlm",
    "spherical2cartesian",
    "sphvec2cart",
    "xlm",
    "ylm",
]

import numpy as np
from scipy.special import factorial, lpmn

try:
    from quadpy import sphere
except ImportError:
    from quadpy import u3 as sphere

from .mesh_calculus import gradient_matrix
from .utils import tri_normals_and_areas


############################################
# Helper functions for sph computation


def cartesian2spherical(p, zaxis_approx=True):
    """
    Maps cartesian coordinates to spherical.

    Parameters
    ----------
    p: Nx3 array
        cartesian coordinates
    zaxis_approx: Boolean (True)
        If True, apply regularization that avoids singularity on z-axis

    Returns
    -------

    sp: Nx3 array
        spherical coordinates [r, theta, phi]

    """

    r = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2 + p[:, 2] ** 2)
    theta = np.arctan2(np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2), p[:, 2])
    phi = np.arctan2(p[:, 1], p[:, 0])

    sp = np.array((r, theta, phi)).T

    if zaxis_approx:
        eps = 1e-6
        theta_mask_0 = np.abs(sp[:, 1]) < eps
        theta_mask_pi = np.abs(sp[:, 1] - np.pi) < eps
        sp[theta_mask_0, 1] += eps
        sp[theta_mask_pi, 1] -= eps

    return sp


def spherical2cartesian(sp):
    """
    Maps spherical coordinates to cartesian.

    Parameters
    ----------
    sp: Nx3 array
        spherical coordinates [r, theta, phi]

    Returns
    -------
    p: Nx3 array
        cartesian croodinates
    """

    X = sp[:, 0] * np.sin(sp[:, 1]) * np.cos(sp[:, 2])
    Y = sp[:, 0] * np.sin(sp[:, 1]) * np.sin(sp[:, 2])
    Z = sp[:, 0] * np.cos(sp[:, 1])

    p = np.array((X, Y, Z)).T
    return p


def Rotmatrix(sp):
    """
    Constructs rotation matrix from cartesian coordinates to spherical.

    Parameters
    ----------
    sp: Nx3 array
        spherical coordinates [r, theta, phi]

    Returns
    -------
    vmat: 3x3 array
        rotation matrix from cartesian to spherical.
    """

    vmat = np.zeros((3, 3))
    vmat[0, 0] = np.sin(sp[1]) * np.cos(sp[2])
    vmat[0, 1] = np.sin(sp[1]) * np.sin(sp[2])
    vmat[0, 2] = np.cos(sp[1])
    vmat[1, 0] = np.cos(sp[1]) * np.cos(sp[2])
    vmat[1, 1] = np.cos(sp[1]) * np.sin(sp[2])
    vmat[1, 2] = -1 * np.sin(sp[1])
    vmat[2, 0] = -1 * np.sin(sp[2])
    vmat[2, 1] = np.cos(sp[2])
    vmat[2, 2] = 0

    return vmat


def cartvec2sph(sp, vec):
    """
    Transforms cartesian vector to spherical coordinates.

    Parameters
    ----------
    sp: Nx3 array
        spherical coordinates [r, theta, phi]
    vec: Nx3 array
        vector in cartesian coordinates [e_x, e_y, e_z]

    Returns
    -------
    svec: Nx3 array
        vector in spherical coordinates [e_r, e_theta, e_phi]

    """

    svec = np.zeros(vec.shape)
    for i in range(sp.shape[0]):
        vmat = Rotmatrix(sp[i])
        svec[i] = vmat @ vec[i]
    return svec


def sphvec2cart(sp, vec):
    """
    Transforms cartesian vector to spherical coordinates.

    Parameters
    ----------
    sp: Nx3 array
        spherical coordinates [r, theta, phi]
    vec: Nx3 array
        vector in spherical coordinates [e_r, e_theta, e_phi]

    Returns
    -------
    svec: Nx3 array
        vector in cartesian coordinates [e_x, e_y, e_z]

    """

    svec = np.zeros(vec.shape)
    for i in range(sp.shape[0]):
        vmat = Rotmatrix(sp[i])
        svec[i] = np.transpose(vmat) @ vec[i]

    return svec


############################################
# Functions for generating Legendre polynomials,
# real spherical harmonics and vector spherical harmonics


def lpmn_em(l, m, x):
    """
    Computes associated Legendre function (Plm) of the first kind of order m and degree l.

    Parameters
    ----------

    l: int
        degree of Plm
    m: int
        order of Plm
    x: Nx1 array
        evaluation points

    Returns
    -------
    lp: Nx1 array
        Plm at `x`

    """

    lp = np.zeros(x.shape)
    for i in range(x.shape[0]):
        a, b = lpmn(m, l, x[i])
        lp[i] = a[np.abs(m), l]
    return lp


def derlpmn_em(l, m, x):
    """
    Computes derivative of associated Legendre function (Plm) of the first kind of order m and degree l
    with respect to the argument x.

    Parameters
    ----------
    l: int
        degree of Plm
    m: int
        order of Plm
    x: Nx1 array
        evaluation points

    Returns
    -------
    derlp: Nx1 array
        dPlm/dx at `x`

    """

    derlp = np.zeros(x.shape)
    for i in range(x.shape[0]):
        a, b = lpmn(m, l, x[i])
        derlp[i] = b[np.abs(m), l]
    return derlp


def xlm(l, m, theta):
    """
    Xlm-function used in the definition of spherical harmonics (Ylm).
    Follows notation of Plattner and Simons (2014);
    see Eqs. 1--3 in https://arxiv.org/pdf/1306.3201.pdf.

    Parameters
    ----------
    l: int
        degree of Xlm
    m: int
        order of Xlm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates

    Returns
    -------
    Xlm: Nx1 array
        Xlm at theta

    """

    xlm = ((2 * l + 1) / (4 * np.pi)) ** 0.5 * (
        factorial(l - m) / factorial(l + m)
    ) ** 0.5
    xlm *= lpmn_em(l, m, np.cos(theta))
    return xlm


def ylm(l, m, theta, phi):
    """
    Real spherical harmonics as defined by Plattner and Simons (2014);
    see Eqs. 1--3 in https://arxiv.org/pdf/1306.3201.pdf.

    Parameters
    ----------
    l: int
        degree of Ylm
    m: int
        order of Ylm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates
    phi: Nx1 array
        evaluation points, phi at spherical coordinates

    Returns
    -------
    Ylm: Nx1 array
        Ylm at (theta,phi)

    """

    if m < 0:
        ylm = np.sqrt(2) * xlm(l, m, theta) * np.cos(m * phi)
    elif m == 0:
        ylm = xlm(l, 0, theta)
    elif m > 0:
        ylm = np.sqrt(2) * xlm(l, m, theta) * np.sin(m * phi)

    return ylm


def derxlm(l, m, theta):
    """
    Derivative of Xlm with respect to theta.

    Parameters
    ----------
    l: int
        degree of Xlm
    m: int
        order of Xlm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates

    Returns
    -------
    derxlm: Nx1 array
        dXlm/dtheta at theta

    """

    derxlm = ((2 * l + 1) / (4 * np.pi)) ** 0.5 * (
        factorial(l - m) / factorial(l + m)
    ) ** 0.5
    derxlm *= derlpmn_em(l, m, np.cos(theta))
    derxlm *= -1 * np.sin(
        theta
    )  # this comes from dXlm(cos(theta))/dtheta = dXlm(cos(theta))/dcos(theta)*(-sin(theta))
    return derxlm


def sinxlm(l, m, theta):
    """
    Computes m/(sin(theta))*Xlm.

    Parameters
    ----------
    l: int
        degree of Xlm
    m: int
        order of Xlm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates

    Returns
    -------
    sinxlm: Nx1 array
        m/(sin(theta))*Xlm at theta

    """

    sinxlm = m / (np.sin(theta)) * xlm(l, m, theta)
    return sinxlm


def dthylm(l, m, theta, phi):
    """
    Derivative of Ylm with respect to theta dYlm/dtheta.

    Parameters
    ----------
    l: int
        degree of Ylm
    m: int
        order of Ylm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates
    phi: Nx1 array
        evaluation points, phi at spherical coordinates

    Returns
    -------
    dthylm: Nx1 array
        dYlm/dtheta at (theta, phi).

    """

    if m < 0:
        dthylm = np.sqrt(2) * derxlm(l, m, theta) * np.cos(m * phi)
    elif m == 0:
        dthylm = derxlm(l, 0, theta)
    elif m > 0:
        dthylm = np.sqrt(2) * derxlm(l, m, theta) * np.sin(m * phi)
    return dthylm


def dphiylm(l, m, theta, phi):
    """
    Derivative of Ylm with respect to phi dYlm/dphi.

    Parameters
    ----------
    l: int
        degree of Ylm
    m: int
        order of Ylm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates
    phi: Nx1 array
        evaluation points, phi at spherical coordinates

    Returns
    -------
    dphiylm: Nx1 array
        dYlm/dphi at (theta, phi).

    """

    if m < 0:
        dphiylm = -np.sqrt(2) * np.sin(m * phi) * sinxlm(l, m, theta)
    elif m == 0:
        dphiylm = 0
    elif m > 0:
        dphiylm = np.sqrt(2) * np.cos(m * phi) * sinxlm(l, m, theta)
    return dphiylm


def Plm(l, m, theta, phi):
    """
    Plm vector function (see Eq. 18 Plattner and Simons (2014)).

    Parameters
    ----------
    l: int
        degree of Plm
    m: int
        order of Plm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates
    phi: Nx1 array
        evaluation points, phi at spherical coordinates

    Returns
    -------
    Plm: Nx3 array
        Plm at (theta, phi).

    """

    Plm = np.zeros((theta.shape[0], 3))
    Plm[:, 0] = ylm(l, m, theta, phi)
    return Plm


def Blm(l, m, theta, phi):
    """
    Blm vector function (see Eq. 19 Plattner and Simons (2014)).

    Parameters
    ----------
    l: int
        degree of Plm
    m: int
        order of Plm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates
    phi: Nx1 array
        evaluation points, phi at spherical coordinates

    Returns
    -------
    Blm: Nx3 array
        Blm at (theta, phi).

    """

    Blm = np.zeros((theta.shape[0], 3))

    Blm[:, 1] = dthylm(l, m, theta, phi)

    Blm[:, 2] = dphiylm(l, m, theta, phi)

    Blm *= 1 / np.sqrt(l * (l + 1))
    return Blm


def Wlm(l, m, theta, phi):
    """
    Vector basis function (Wlm) for r**l component of the magnetic field.
    Normalization <Wlm,Wkn> = delta_lk,mn.

    Parameters
    ----------
    l: int
        degree of Wlm
    m: int
        order of Wlm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates
    phi: Nx1 array
        evaluation points, phi at spherical coordinates

    Returns
    -------
    Wlm: Nx3 array
        Wlm at (theta, phi).

    """

    Wlm = l * Plm(l, m, theta, phi) + np.sqrt(l * (l + 1)) * Blm(l, m, theta, phi)
    Wlm *= 1 / np.sqrt(2 * l ** 2 + l)
    return Wlm


def Vlm(l, m, theta, phi):
    """
    Vector basis function (Vlm) for r**(-l) component of the magnetic field.
    Normalization <Vlm,Vkn> = delta_lk,mn.

    Parameters
    ----------
    l: int
        degree of Vlm
    m: int
        order of Vlm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates
    phi: Nx1 array
        evaluation points, phi at spherical coordinates

    Returns
    -------
    Vlm: Nx3 array
        Vlm at (theta, phi).

    """

    Vlm = -1 * (l + 1) * Plm(l, m, theta, phi) + np.sqrt(l * (l + 1)) * Blm(
        l, m, theta, phi
    )
    Vlm *= 1 / np.sqrt((l + 1) * (2 * l + 1))
    return Vlm


def Xlm(l, m, theta, phi):
    """
    Vector spherical harmonics basis function (Xlm).
    Normalization <Xlm,Xkn> = delta_lk,mn.

    Parameters
    ----------
    l: int
        degree of Xlm
    m: int
        order of Xlm
    theta: Nx1 array
        evaluation points, theta at spherical coordinates
    phi: Nx1 array
        evaluation points, phi at spherical coordinates

    Returns
    -------
    Xlm: Nx3 array
        Xlm at (theta, phi).

    """

    temp = Blm(l, m, theta, phi)
    r = np.zeros(temp.shape)
    r[:, 0] = 1
    Xlm = -1 * np.cross(r, temp, axis=1)
    return Xlm


############################################
# Potential and field from the spherical harmonics


def basis_potentials(p, lmax, normalization="default", R=1):
    """
    Computes inner/outer basis functions for magnetic scalar potential.
    Ignores the 'DC' component l=0.

    Parameters
    ----------
    p: Nx3 array
        coordinates in which the potential is computed
    lmax: int
        maximum degree l which is used in computing
    normalization: string
        which normalization scheme to use
        "default": Integration over unit sphere <Y_lm,Y_l'm'> = delta_ll'mm' (default)
        "energy": each term in inner/outer basis function with respect to 
        sphere with radius R is normalized to unit energy
    R: float
        Sphere radius that separates inner/outer components            
                  
    Returns
    -------
    pot1: Nxlmax*(lmax+2) array
        inner potential (alpha_lm) basis functions at p
    pot2: Nxlmax*(lmax+2) array
        outer potential (beta_lm) basis functions at p     

    """
    mu0 = 1e-7 * 4 * np.pi
    L = lmax * (lmax + 2)

    pot1 = np.zeros((p.shape[0], L))
    pot2 = np.zeros((p.shape[0], L))

    sp = cartesian2spherical(p)

    lind = 0
    for l in range(1, lmax + 1):
        for m in range(-1 * l, l + 1):
            _ylm = ylm(l, m, sp[:, 1], sp[:, 2])

            pot1[:, lind] = sp[:, 0] ** (-l - 1) * _ylm
            pot2[:, lind] = sp[:, 0] ** l * _ylm

            if normalization == "energy":
                pot1[:, lind] *= np.sqrt(R ** (2 * l + 1) / ((l + 1) * mu0))
                pot2[:, lind] *= 1 / np.sqrt(R ** (2 * l + 1) * l * mu0)

            lind += 1

    return pot1, pot2


def potential(p, acoeffs, bcoeffs, lmax, normalization="default", R=1):
    """
    Computes magnetic scalar potential from the sph coefficients.
    Ignores the 'DC' component l=0.

    Parameters
    ----------
    p: Nx3 array
        coordinates in which the potential is computed
    acoeffs: lmax*(lmax+2)x1 array
        spectral coefficients of r**l terms
    bcoeffs: lmax*(lmax+2)x1 array
        spectral coefficients of r**(-l) terms
    lmax: int
        maximum degree l which is used in computing
    normalization: string
        which normalization scheme to use
        "default": Integration over unit sphere <Y_lm,Y_l'm'> = delta_ll'mm' (default)
        "energy": each term in inner/outer basis function with respect to 
        sphere with radius R is normalized to unit energy
    R: float
        Sphere radius that separates inner/outer components   
                  
    Returns
    -------
    pot1: Nx array
        magnetic scalar potential at p

    """

    basis = basis_potentials(p, lmax, normalization, R)

    return basis[0] @ acoeffs + basis[1] @ bcoeffs


def field(p, acoeffs, bcoeffs, lmax, normalization="default", R=1):
    """
      Computes magnetic field at some point from the sph coefficients.
    Ignores the 'DC' component l=0.

    Parameters
    ----------
    p: Nx3 array
        coordinates in which the potential is computed
    acoeffs: lmax*(lmax+2)x1 array
        spectral coefficients of r**l terms
    bcoeffs: lmax*(lmax+2)x1 array
        spectral coefficients of r**(-l) terms
    lmax: int
        maximum degree l which is used in computing  
    normalization: string
        "default": the fields correspond to normalized magnetic potential.
        "unit": the fields are normalized w.r.t integration over the unit sphere.
        "energy": each term in inner/outer basis function with respect to 
        sphere with radius R is normalized to unit energy 
    R: float
        Sphere radius that separates inner/outer components   
                
    Returns
    -------
    B: Nx3 array
        Magnetic field produced by the sph components
    """

    basis = basis_fields(p, lmax, normalization=normalization, R=R)

    return basis[0] @ acoeffs + basis[1] @ bcoeffs


def basis_fields(p, lmax, normalization="default", R=1):
    """
    Computes magnetic field of each sph coefficient.
    Ignores the 'DC' component l=0.

    Parameters
    ----------
    p: Nx3 array
        coordinates in which the field is computed
    lmax: int
        maximum degree l used in the computation
    normalization: string
        "default": the fields correspond to magnetic normalized potential.
        "unit": the fields are normalized w.r.t integration over the unit sphere.
        "energy": each term in inner/outer basis function with respect to 
            sphere with radius R is normalized to unit energy 
    R: float
        Sphere radius that separates inner/outer components   
                
    Returns
    -------
    B1: N x 3 x N_lmax array
        magnetic field at p for each alpha_lm
    B2: N x 3 x N_lmax array
        magnetic field at p for each beta_lm

    """
    mu0 = 1e-7 * 4 * np.pi
    L = lmax * (lmax + 2)  # Fixed
    B1 = np.zeros((L, p.shape[0], p.shape[1]))
    B2 = np.zeros((L, p.shape[0], p.shape[1]))

    sp = cartesian2spherical(p)

    idx = 0
    for l in range(1, lmax + 1):
        for m in range(-1 * l, l + 1):
            _Wlm = Wlm(l, m, sp[:, 1], sp[:, 2])
            _Wlm[:, 0] *= sp[:, 0] ** (l - 1)
            _Wlm[:, 1] *= sp[:, 0] ** (l - 1)
            _Wlm[:, 2] *= sp[:, 0] ** (l - 1)

            if normalization == "default":
                _Wlm *= np.sqrt(2 * l ** 2 + l) * mu0
            if normalization == "energy":
                _Wlm *= np.sqrt(2 * l ** 2 + l) * mu0
                _Wlm *= 1 / np.sqrt(R ** (2 * l + 1) * l * mu0)

            _Wlm = sphvec2cart(sp, _Wlm)
            B2[idx] = -1 * _Wlm  # r**l functions

            _Vlm = Vlm(l, m, sp[:, 1], sp[:, 2])
            if normalization == "default":
                _Vlm *= np.sqrt((2 * l + 1) * (l + 1)) * mu0
            if normalization == "energy":
                _Vlm *= np.sqrt((2 * l + 1) * (l + 1)) * mu0
                _Vlm *= np.sqrt(R ** (2 * l + 1) / ((l + 1) * mu0))

            _Vlm[:, 0] *= sp[:, 0] ** (-l - 2)
            _Vlm[:, 1] *= sp[:, 0] ** (-l - 2)
            _Vlm[:, 2] *= sp[:, 0] ** (-l - 2)
            _Vlm = sphvec2cart(sp, _Vlm)
            B1[idx] = -1 * _Vlm  # 1/r**l functions

            idx += 1

    # FIX, should be handled earlier maybe
    B1[np.isinf(B1)] = 0
    B2[np.isinf(B2)] = 0

    return np.moveaxis(B1, 0, 2), np.moveaxis(B2, 0, 2)


def compute_sphcoeffs_mesh(mesh, lmax, normalization="default", R=1):
    """
    Computes multipole moment (spherical harmonics coefficient) transformation
    from the mesh.

    Parameters
    ----------
    mesh: mesh object - the surface mesh
    lmax: int
        maximum degree l of the fit
    normalization: str
        'default' (Ylm**2 integrated over solid angle to 1)
        'energy' (field energy of basis fields normalized to 1 in R-ball)
    R: float                           
        radius in the energy normalization
    Returns
    -------
    Calm: (lmax*(lmax+2)xNvertices array
          transformation from the mesh to alm coefficients (r**l-terms)
    Cblm: (lmax*(lmax+2)xNvertices array
          transformation from the mesh to blm coefficients (r**(-l)-terms)

    """
    mu0 = 1e-7 * 4 * np.pi

    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)
    tri_normals, tri_areas = tri_normals_and_areas(mesh.vertices, mesh.faces)

    centers = np.mean(mesh.vertices[mesh.faces], axis=1)
    centers_sp = cartesian2spherical(centers)

    Calm = np.zeros((lmax * (lmax + 2), mesh.vertices.shape[0]))
    Cblm = np.zeros((lmax * (lmax + 2), mesh.vertices.shape[0]))

    idx = 0
    for l in range(1, lmax + 1):
        for m in range(-1 * l, l + 1):
            _Xlm = sphvec2cart(centers_sp, Xlm(l, m, *centers_sp[:, 1:].T))
            # Multiply by the normalization factor, and tri areas
            _Xlm *= np.sqrt(l * (l + 1)) * tri_areas[:, None]
            # alm ja blm have different r-dependencies
            alm_terms = _Xlm.T * (centers_sp[:, 0] ** (l))
            blm_terms = _Xlm.T * (centers_sp[:, 0] ** (-l - 1))

            # Combine rotated gradient (operates on stream function -> j)
            integral_alm = alm_terms[0] @ Gx + alm_terms[1] @ Gy + alm_terms[2] @ Gz
            integral_blm = blm_terms[0] @ Gx + blm_terms[1] @ Gy + blm_terms[2] @ Gz
            # Default normalization
            Cblm[idx] = -1 / ((2 * l + 1) * l) * integral_blm
            Calm[idx] = 1 / ((2 * l + 1) * (l + 1)) * integral_alm
            if normalization == "energy":
                Cblm[idx] *= np.sqrt(R ** (2 * l + 1) * l * mu0)
                Calm[idx] *= np.sqrt(((l + 1) * mu0) / (R ** (2 * l + 1)))

            idx += 1
        print("l = %d computed" % (l))

    return Calm, Cblm


###################################
# sph class


class SphBasis:
    """
    Constructs an object that describes unit sphere. Can be used to compute
    inner products and function sph spectra on the sphere.

    TODO: mu0 might be missing!!!

    Properties
    ----------
    sp: Nx3 array
        spherical coordinates of the points [r, theta, phi]
    p: Nx3 array
        cartesian coordinates
    qp: Mx3 array
        cartesian quadrature points
    sqp: Mx3 array
        spherical quadrature points
    """

    def __init__(self, Np):
        """
        Initialises the sphbasis object.

        Parameters
        ----------
        Np: int
            Mumber of points along theta and phi in spherical coordinates

        Returns
        -------
        self: sphbasis object

        """

        theta = np.linspace(0.01, np.pi - 0.01, Np)
        phi = np.linspace(0, 2 * np.pi, Np)
        phi, theta = np.meshgrid(phi, theta)
        phi = phi.flatten()
        theta = theta.flatten()

        self.sp = np.zeros((theta.shape[0], 3))
        self.sp[:, 0] = 1
        self.sp[:, 1] = theta
        self.sp[:, 2] = phi

        self.p = spherical2cartesian(self.sp)
        self.Np = Np

        self.initqp()

    def initqp(self):
        """
        Initialises quadrature points on the sphere.

        Default points are McLaren(10) so that we avoid singularities.
        """

        self.qp = sphere.mclaren_10()
        sp = cartesian2spherical(self.qp.points)
        self.sqp = sp

    def innerproduct(self, fun1, fun2):
        """
        Inner product of vector functions fun1 and fun2.
        Defined as integration over a surface of unit sphere <C,D> = int C dot D dOmega.
        Quadrature rule defined in qp is used.

        Parameters
        ----------
        fun1: Nx3 array
            vector function 1
        fun2: Nx3 array
            vector function 2

        Returns
        -------
        dotp: int
            inner product of fun1 and fun2

        """

        dotp = np.sum(self.qp.weights * np.sum(fun1 * fun2, axis=1)) * 4 * np.pi
        return dotp

    def avsphspectra(self, fun, lmax):
        """
        Calculate the l,m-spectra (over r**l-terms) of vector function defined in quadrature points
        using the inner product.

        Parameters
        ----------
        fun: Nx3 array
            vector function computed at quadrature points self.sqp
        lmax: int
            maximum degree l for which the spectra is computed

        Returns
        -------
        coeffs: lmax*(lmax+2)x1 arrays
            spectral coefficients

        """

        coeffs = []

        for l in range(1, lmax + 1):
            for m in range(-1 * l, l + 1):
                _Wlm = Wlm(l, m, self.sqp[:, 1], self.sqp[:, 2])
                ctemp = self.innerproduct(fun, _Wlm)
                ctemp /= self.sqp[0, 0] ** (l - 1) * np.sqrt(
                    2 * l ** 2 + l
                )  # we use this normalization
                #                ctemp /= (self.sqp[0,0]**(l-1))
                coeffs.append(ctemp)

        coeffs = np.array(coeffs)
        return coeffs

    def bvsphspectra(self, fun, lmax):
        """
        Calculate the l,m-spectra (over r**(-l)-terms) of vector function defined in quadrature points
        using the inner product.

        Parameters
        ----------
        fun: Nx3 array
            vector function computed at quadrature points self.sqp
        lmax: int
            maximum degree l for which the spectra is computed

        Returns
        -------
        coeffs: lmax*(lmax+2)x1 arrays
            spectral coefficients

        """

        coeffs = []

        for l in range(1, lmax + 1):
            for m in range(-1 * l, l + 1):
                _Vlm = Vlm(l, m, self.sqp[:, 1], self.sqp[:, 2])
                ctemp = self.innerproduct(fun, _Vlm)
                ctemp /= self.sqp[0, 0] ** (l - 1) * np.sqrt(
                    (l + 1) * (2 * l + 1)
                )  # we use this normalization
                #                ctemp /= (self.sqp[0,0]**(l-1))
                coeffs.append(ctemp)

        coeffs = np.array(coeffs)
        return coeffs


###################################
# Functions for plotting spherical harmonics


def plotYlms(sph, lmax, polar=False):
    """
    Plots real spherical harmonics up to lmax.
    Inspired by https://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html.

    Parameters
    ----------
    sph: spherical harmonics analysis object
    lmax: int

        maximum degree l
    polar: boolean
        plot polar representation?
    """
    from mayavi import mlab

    theta = np.reshape(sph.sp[:, 1], (sph.Np, sph.Np))
    phi = np.reshape(sph.sp[:, 2], (sph.Np, sph.Np))
    r = 0.4
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    if polar:
        for l in range(1, lmax + 1):
            for m in range(l):
                _ylm = ylm(l, m, theta.flatten(), phi.flatten())
                _ylm = np.reshape(_ylm, (sph.Np, sph.Np))

                mlab.mesh(x - m, y - l, z, scalars=_ylm, colormap="bwr")
                _ylm /= _ylm.max()
                mlab.mesh(
                    _ylm * x - m,
                    _ylm * y - l,
                    _ylm * z + 1.3,
                    scalars=np.abs(_ylm),
                    colormap="Spectral",
                )

        mlab.view(90, 70, 6.2, (-1.3, -2.9, 0.25))
    else:
        for l in range(0, lmax + 1):
            for m in range(-l, l + 1):
                _ylm = ylm(l, m, theta.flatten(), phi.flatten())
                _ylm = np.reshape(_ylm, (sph.Np, sph.Np))

                mlab.mesh(x - m, y - l, z, scalars=_ylm, colormap="bwr")

        mlab.view(0, 180)


def plotYlm(sph, l, m):
    """
    Plots real spherical harmonics of order m and degree l.
    Inspired by https://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html

    Parameters
    ----------
    sph: spherical harmonics analysis object
    l: int
        degree l
    m: int
        order m

    """
    from mayavi import mlab

    theta = np.reshape(sph.sp[:, 1], (sph.Np, sph.Np))
    phi = np.reshape(sph.sp[:, 2], (sph.Np, sph.Np))
    r = 0.6
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    _ylm = ylm(l, m, theta.flatten(), phi.flatten())
    _ylm = np.reshape(_ylm, (sph.Np, sph.Np))

    mlab.mesh(x - m, y - l, z, scalars=_ylm, colormap="bwr")

    _ylm /= _ylm.max()
    mlab.mesh(
        _ylm * x - m,
        _ylm * y - l,
        _ylm * z + 1.3,
        scalars=np.abs(_ylm),
        colormap="Spectral",
    )


def plotWlm(sph, l, m):
    """
    Plots magnetic field basis function 'Wlm' (r**l) over a sphere.

    Parameters
    ----------
    sph: spherical harmonics analysis object
    l: int
        degree l
    m: int
        order m

    Returns
    -------
    obj: mayavi object

    """
    from mayavi import mlab

    _Wlm = Wlm(l, m, sph.sp[:, 1], sph.sp[:, 2])
    _Wlm = sphvec2cart(sph.sp, _Wlm)
    obj = mlab.quiver3d(
        sph.p[:, 0], sph.p[:, 1], sph.p[:, 2], _Wlm[:, 0], _Wlm[:, 1], _Wlm[:, 2]
    )
    obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
    return obj


def plotBWlm_volume(sph, l, m, lim, Np, offset):
    """
    Plots magnetic field basis function 'Wlm' (r**l) over a 3D volume.

    Parameters
    ----------
    sph: spherical harmonics analysis object
    l: int
        degree l
    m: int
        order m
    lim: float
        limits for coordinates, e.g., xmin = -lim, xmax = lim
    Np: int
        number of points along different coordinates
    offset: 1x3 array
        offset of the volume in which Wlm is plotted

    Returns
    -------
        obj: mayavi object

    """

    from mayavi import mlab

    x, y, z = np.meshgrid(
        np.linspace(-lim + offset[0], lim + offset[0], Np),
        np.linspace(-lim + offset[1], lim + offset[1], Np),
        np.linspace(-lim + offset[2], lim + offset[2], Np),
    )

    p = np.array((x.flatten(), y.flatten(), z.flatten())).T
    sp = cartesian2spherical(p)

    _Wlm = Wlm(l, m, sp[:, 1], sp[:, 2])
    _Wlm *= np.sqrt(2 * l ** 2 + l)
    _Wlm[:, 0] *= sp[:, 0] ** (l - 1)
    _Wlm[:, 1] *= sp[:, 0] ** (l - 1)
    _Wlm[:, 2] *= sp[:, 0] ** (l - 1)

    _Wlm = sphvec2cart(sp, _Wlm)
    obj = mlab.quiver3d(p[:, 0], p[:, 1], p[:, 2], _Wlm[:, 0], _Wlm[:, 1], _Wlm[:, 2])
    obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
    return obj


def plotVlm(sph, l, m):
    """
    Plots magnetic field basis function 'Vlm' (r**(-l)) over a sphere.

    Parameters
    ----------
    sph: spherical harmonics analysis object
    l: int
        degree l
    m: int
        order m

    Returns
    -------
    obj: mayavi object

    """
    from mayavi import mlab

    _Vlm = Vlm(l, m, sph.sp[:, 1], sph.sp[:, 2])
    _Vlm = sphvec2cart(sph.sp, _Vlm)
    obj = mlab.quiver3d(
        sph.p[:, 0], sph.p[:, 1], sph.p[:, 2], _Vlm[:, 0], _Vlm[:, 1], _Vlm[:, 2]
    )
    obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
    return obj


def plotBVlm_volume(sph, l, m, lim, Np, offset):
    """
    Plots magnetic field basis function 'Vlm' (r**(-l)) over a 3D volume.

    Parameters
    ----------
    sph: spherical harmonics analysis object
    l: int
        degree l
    m: int
        order m
    lim: float
        limits for coordinates, e.g., xmin = -lim, xmax = lim
    Np: int
        number of points along different coordinates
    offset: 1x3 array
        offset of the volume in which Vlm is plotted

    Returns
    -------
    obj: mayavi object

    """
    from mayavi import mlab

    x, y, z = np.meshgrid(
        np.linspace(-lim + offset[0], lim + offset[0], Np),
        np.linspace(-lim + offset[1], lim + offset[1], Np),
        np.linspace(-lim + offset[2], lim + offset[2], Np),
    )

    p = np.array((x.flatten(), y.flatten(), z.flatten())).T
    sp = cartesian2spherical(p)

    _Vlm = Vlm(l, m, sp[:, 1], sp[:, 2])
    _Vlm *= np.sqrt((l + 1) * (2 * l + 1))

    _Vlm[:, 0] *= sp[:, 0] ** (-1 * (l + 2))
    _Vlm[:, 1] *= sp[:, 0] ** (-1 * (l + 2))
    _Vlm[:, 2] *= sp[:, 0] ** (-1 * (l + 2))

    _Vlm = sphvec2cart(sp, _Vlm)
    obj = mlab.quiver3d(p[:, 0], p[:, 1], p[:, 2], _Vlm[:, 0], _Vlm[:, 1], _Vlm[:, 2])
    obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
    return obj


def plotXlm(sph, l, m):
    """
    Plots vector spherical harmonic basis function 'Xlm' over a sphere.

    Parameters
    ----------
    sph: spherical harmonics analysis object
    l: int
        degree l
    m: int
        order m

    Returns
    -------
    obj: mayavi object

    """
    from mayavi import mlab

    _Xlm = Xlm(l, m, sph.sp[:, 1], sph.sp[:, 2])
    _Xlm = sphvec2cart(sph.sp, _Xlm)
    obj = mlab.quiver3d(
        sph.p[:, 0], sph.p[:, 1], sph.p[:, 2], _Xlm[:, 0], _Xlm[:, 1], _Xlm[:, 2]
    )
    obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
    return obj
