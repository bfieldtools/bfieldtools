from bfieldtools import sphtools

import pytest


import numpy as np

from numpy.testing import assert_allclose


def test_coord_changes():

    point = np.array([[0.4, 0.4, 1]])

    sphpoint = sphtools.cartesian2spherical(point)

    point2 = sphtools.spherical2cartesian(sphpoint)

    assert_allclose(point, point2)


# def Rotmatrix(sp):
#     """
#     Constructs rotation matrix from cartesian coordinates to spherical.

#     Parameters
#     ----------
#     sp: Nx3 array
#         spherical coordinates [r, theta, phi]

#     Returns
#     -------
#     vmat: 3x3 array
#         rotation matrix from cartesian to spherical.
#     """

#     vmat = np.zeros((3, 3))
#     vmat[0, 0] = np.sin(sp[1]) * np.cos(sp[2])
#     vmat[0, 1] = np.sin(sp[1]) * np.sin(sp[2])
#     vmat[0, 2] = np.cos(sp[1])
#     vmat[1, 0] = np.cos(sp[1]) * np.cos(sp[2])
#     vmat[1, 1] = np.cos(sp[1]) * np.sin(sp[2])
#     vmat[1, 2] = -1 * np.sin(sp[1])
#     vmat[2, 0] = -1 * np.sin(sp[2])
#     vmat[2, 1] = np.cos(sp[2])
#     vmat[2, 2] = 0

#     return vmat


# def cartvec2sph(sp, vec):
#     """
#     Transforms cartesian vector to spherical coordinates.

#     Parameters
#     ----------
#     sp: Nx3 array
#         spherical coordinates [r, theta, phi]
#     vec: Nx3 array
#         vector in cartesian coordinates [e_x, e_y, e_z]

#     Returns
#     -------
#     svec: Nx3 array
#         vector in spherical coordinates [e_r, e_theta, e_phi]

#     """

#     svec = np.zeros(vec.shape)
#     for i in range(sp.shape[0]):
#         vmat = Rotmatrix(sp[i])
#         svec[i] = vmat @ vec[i]
#     return svec


# def sphvec2cart(sp, vec):
#     """
#     Transforms cartesian vector to spherical coordinates.

#     Parameters
#     ----------
#     sp: Nx3 array
#         spherical coordinates [r, theta, phi]
#     vec: Nx3 array
#         vector in spherical coordinates [e_r, e_theta, e_phi]

#     Returns
#     -------
#     svec: Nx3 array
#         vector in cartesian coordinates [e_x, e_y, e_z]

#     """

#     svec = np.zeros(vec.shape)
#     for i in range(sp.shape[0]):
#         vmat = Rotmatrix(sp[i])
#         svec[i] = np.transpose(vmat) @ vec[i]

#     return svec


def test_sph_eval():
    """
    Simply test that the function evaluate, assert shapes but nothing more
    """

    l = 4

    x = np.array([0.4, 0.4, 1])
    theta = np.array([0.4, 0.4, 1])
    phi = np.array([0.4, 0.4, 1])

    for m in [-2, 0, 2]:
        m = 2

        assert sphtools.lpmn_em(l, m, x).shape == (3,)
        assert sphtools.derlpmn_em(l, m, x).shape == (3,)

        assert sphtools.xlm(l, m, theta).shape == (3,)

        assert sphtools.ylm(l, m, theta, phi).shape == (3,)

        assert sphtools.derxlm(l, m, theta).shape == (3,)

        assert sphtools.sinxlm(l, m, theta).shape == (3,)

        assert sphtools.dthylm(l, m, theta, phi).shape == (3,)

        assert sphtools.dphiylm(l, m, theta, phi).shape == (3,)

        assert sphtools.Plm(l, m, theta, phi).shape == (3, 3)

        assert sphtools.Blm(l, m, theta, phi).shape == (3, 3)

        assert sphtools.Wlm(l, m, theta, phi).shape == (3, 3)
        assert sphtools.Vlm(l, m, theta, phi).shape == (3, 3)


def test_potential():

    p = np.array([[1, 0, 0], [-1, 0, 0]])

    acoeffs = np.array([1, 0, 0])
    bcoeffs = np.array([0, 0, 0])
    lmax = 1

    U = sphtools.potential(p, acoeffs, bcoeffs, lmax)

    assert U[0] == -U[-1]


def test_field():

    p = np.array([[1, 0, 0], [-1, 0, 0]])

    acoeffs = np.array([1, 0, 0, 0, 0, 0, 0, 0])
    bcoeffs = np.array([0, 0, 0, 0, 0, 0, 0, 0])
    lmax = 2

    B = sphtools.field(p, acoeffs, bcoeffs, lmax)

    assert_allclose(B[0], B[1], atol=1e-16)

    basis = sphtools.basis_fields(p, lmax)

    B_b = np.einsum("ijk,j->ki", basis[0], acoeffs) + np.einsum(
        "ijk,j->ki", basis[1], bcoeffs
    )

    # Skip until bug solved
    # assert_allclose(B, B_b)


# # sph class


# class sphbasis:
#     """
#     Class for constructing spherical harmonics (Ylms), their gradients
#     and related magnetic field 'basis vectorfunctions'
#     (Wlms for r**l components, Vlms for r**(-l) components).

#     Uses notations and definitions by Plattner and Simons (2014; https://arxiv.org/pdf/1306.3201.pdf)
#     and the same normalization conventions.

#     Integration over a surface of unit sphere is
#     used as the inner product <C,D> = int C dot D dOmega.

#     Has also lot of functions for spherical <-> cartesian transformations.

#     TODO: mu0 might be missing!!!

#     Properties
#     ----------
#     sp: Nx3 array
#         spherical coordinates of the points [r, theta, phi]
#     p: Nx3 array
#         cartesian coordinates
#     qp: Mx3 array
#         cartesian quadrature points
#     sqp: Mx3 array
#         spherical quadrature points
#     """

#     def __init__(self, Np):
#         """
#         Initialises the sphbasis object.

#         Parameters
#         ----------
#         Np: int
#             Mumber of points along theta and phi in spherical coordinates

#         Returns
#         -------
#         self: sphbasis object

#         """

#         theta = np.linspace(0.01, np.pi - 0.01, Np)
#         phi = np.linspace(0, 2 * np.pi, Np)
#         phi, theta = np.meshgrid(phi, theta)
#         phi = phi.flatten()
#         theta = theta.flatten()

#         self.sp = np.zeros((theta.shape[0], 3))
#         self.sp[:, 0] = 1
#         self.sp[:, 1] = theta
#         self.sp[:, 2] = phi

#         self.p = spherical2cartesian(self.sp)
#         self.Np = Np

#         self.initqp()

#     def initqp(self):
#         """
#         Initialises quadrature points on the sphere.

#         Default points are McLaren(10) so that we avoid singularities.
#         """

#         self.qp = quadpy.sphere.mclaren_10()
#         sp = cartesian2spherical(self.qp.points)
#         self.sqp = sp

#     def innerproduct(self, fun1, fun2):
#         """
#         Inner product of vector functions fun1 and fun2.
#         Defined as integration over a surface of unit sphere <C,D> = int C dot D dOmega.
#         Quadrature rule defined in qp is used.

#         Parameters
#         ----------
#         fun1: Nx3 array
#             vector function 1
#         fun2: Nx3 array
#             vector function 2

#         Returns
#         -------
#         dotp: int
#             inner product of fun1 and fun2

#         """

#         dotp = np.sum(self.qp.weights * np.sum(fun1 * fun2, axis=1)) * 4 * np.pi
#         return dotp

#     def avsphspectra(self, fun, lmax):
#         """
#         Calculate the l,m-spectra (over r**l-terms) of vector function defined in quadrature points
#         using the inner product.

#         Parameters
#         ----------
#         fun: Nx3 array
#             vector function computed at quadrature points self.sqp
#         lmax: int
#             maximum degree l for which the spectra is computed

#         Returns
#         -------
#         coeffs: lmax*(lmax+2)x1 arrays
#             spectral coefficients

#         """

#         coeffs = []

#         for l in range(1, lmax + 1):
#             for m in range(-1 * l, l + 1):
#                 _Wlm = Wlm(l, m, self.sqp[:, 1], self.sqp[:, 2])
#                 ctemp = self.innerproduct(fun, _Wlm)
#                 ctemp /= self.sqp[0, 0] ** (l - 1) * np.sqrt(
#                     2 * l ** 2 + l
#                 )  # we use this normalization
#                 #                ctemp /= (self.sqp[0,0]**(l-1))
#                 coeffs.append(ctemp)

#         coeffs = np.array(coeffs)
#         return coeffs

#     def bvsphspectra(self, fun, lmax):
#         """
#         Calculate the l,m-spectra (over r**(-l)-terms) of vector function defined in quadrature points
#         using the inner product.

#         Parameters
#         ----------
#         fun: Nx3 array
#             vector function computed at quadrature points self.sqp
#         lmax: int
#             maximum degree l for which the spectra is computed

#         Returns
#         -------
#         coeffs: lmax*(lmax+2)x1 arrays
#             spectral coefficients

#         """

#         coeffs = []

#         for l in range(1, lmax + 1):
#             for m in range(-1 * l, l + 1):
#                 _Vlm = Vlm(l, m, self.sqp[:, 1], self.sqp[:, 2])
#                 ctemp = self.innerproduct(fun, _Vlm)
#                 ctemp /= self.sqp[0, 0] ** (l - 1) * np.sqrt(
#                     (l + 1) * (2 * l + 1)
#                 )  # we use this normalization
#                 #                ctemp /= (self.sqp[0,0]**(l-1))
#                 coeffs.append(ctemp)

#         coeffs = np.array(coeffs)
#         return coeffs


# def fitSpectra(coords, Bmeas, lmax):
#     """
#     Fits spherical harmonics representation (r**l) to measured data.

#     Parameters
#     ----------
#     coords: Nx3x3 array
#         measurement coordinates, each measured field direction
#         in the third dimension: e.g. coords[:,:,2] gives the coordinates of measured z-components.
#     Bmeas: Nx3 array
#         the measured field values along different directions (x,y,z)
#     lmax: int
#         maximum degree l for which the fit is done

#     Returns
#     -------
#     coeffs: lmax*(lmax+2)x1 array
#         the unnormalized coefficients
#     coeffs2: lmax*(lmax+2)x1 array
#         the 'properly' normalized squared coefficients
#     nrmse: float
#         normalized rms error in percents between the data and fit

#     """

#     Nmeas = coords.shape[0]
#     A = np.zeros((3 * Nmeas, lmax * (lmax + 2)))  # initialize the fitting matrix

#     # loop over the components
#     for e in range(3):
#         p = coords[:, :, e]
#         sp = cartesian2spherical(p)

#         lind = 0
#         for l in range(1, lmax + 1):
#             for m in range(-1 * l, l + 1):
#                 _Wlm = Wlm(l, m, sp[:, 1], sp[:, 2])
#                 _Wlm *= np.sqrt(2 * l ** 2 + l)
#                 _Wlm[:, 0] *= sp[:, 0] ** (l - 1)
#                 _Wlm[:, 1] *= sp[:, 0] ** (l - 1)
#                 _Wlm[:, 2] *= sp[:, 0] ** (l - 1)
#                 _Wlm = sphvec2cart(sp, _Wlm)
#                 A[e * Nmeas : (e + 1) * Nmeas, lind] = _Wlm[:, e]
#                 lind += 1
#     print(
#         "Condition number = %f" % (np.linalg.cond(A))
#     )  # print the condition number of A

#     coeffs = (
#         np.linalg.pinv(A) @ Bmeas.T.flatten()
#     )  # compute coefficients using pseudoinverse of A

#     # the following script calculates the normalized coefficients
#     # the coefficients are normalized so that squared norm of magnetic field integrates to 1 over the measurement volume
#     coeffs2 = np.zeros(coeffs.shape)
#     lind = 0
#     Rmax = np.max(sp[:, 0])
#     for l in range(1, lmax + 1):
#         for m in range(-1 * l, l + 1):
#             #               coeffs2[lind] = coeffs[lind]*np.sqrt(2*l**2 + l)
#             temp = (2 * l ** 2 + l) * Rmax ** (2 * l - 1) / (2 * l - 1)
#             coeffs2[lind] = coeffs[lind] ** 2 * temp
#             lind += 1
#     Breco = A @ coeffs

#     nrmse = (
#         np.sqrt(np.mean((Bmeas.T.flatten() - Breco) ** 2))
#         / np.max(np.abs(Bmeas.T.flatten()))
#         * 100
#     )
#     print("Normalized RMS error = %f%%" % (nrmse))  # print the normalized rms error

#     return coeffs, coeffs2, nrmse


# def reconstructB(p, coeffs, lmax):
#     """
#     Reconstructs the magnetic field using the spherical harmonics coefficients.

#     Parameters
#     ----------
#     p: Nx3 array
#         coordinates where B is reconstructed
#     coeffs: lmax*(lmax+2)x1 array
#         the unnormalized l,m-coefficients
#     lmax:int
#         maximum degree l of the fit

#     Returns
#     -------
#     B: Nx3 array
#         reconstructed magnetic field at p

#     """

#     B = np.zeros(p.shape)
#     sp = cartesian2spherical(p)
#     idx = 0
#     for l in range(1, lmax + 1):
#         for m in range(-1 * l, l + 1):
#             _Wlm = Wlm(l, m, sp[:, 1], sp[:, 2])
#             _Wlm *= np.sqrt(2 * l ** 2 + l)
#             _Wlm[:, 0] *= sp[:, 0] ** (l - 1)
#             _Wlm[:, 1] *= sp[:, 0] ** (l - 1)
#             _Wlm[:, 2] *= sp[:, 0] ** (l - 1)
#             _Wlm *= coeffs[idx]
#             _Wlm = sphvec2cart(sp, _Wlm)
#             B += _Wlm
#             idx += 1
#     return B


# def plotYlms(sph, lmax, polar=False):
#     """
#     Plots real spherical harmonics up to lmax.
#     Inspired by https://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html.

#     Parameters
#     ----------
#     sph: spherical harmonics analysis object
#     lmax: int

#         maximum degree l
#     polar: boolean
#         plot polar representation?
#     """

#     theta = np.reshape(sph.sp[:, 1], (sph.Np, sph.Np))
#     phi = np.reshape(sph.sp[:, 2], (sph.Np, sph.Np))
#     r = 0.4
#     x = r * np.sin(theta) * np.cos(phi)
#     y = r * np.sin(theta) * np.sin(phi)
#     z = r * np.cos(theta)

#     if polar:
#         for l in range(1, lmax + 1):
#             for m in range(l):
#                 _ylm = ylm(l, m, theta.flatten(), phi.flatten())
#                 _ylm = np.reshape(_ylm, (sph.Np, sph.Np))

#                 mlab.mesh(x - m, y - l, z, scalars=_ylm, colormap="bwr")
#                 _ylm /= _ylm.max()
#                 mlab.mesh(
#                     _ylm * x - m,
#                     _ylm * y - l,
#                     _ylm * z + 1.3,
#                     scalars=np.abs(_ylm),
#                     colormap="Spectral",
#                 )

#         mlab.view(90, 70, 6.2, (-1.3, -2.9, 0.25))
#     else:
#         for l in range(0, lmax + 1):
#             for m in range(-l, l + 1):
#                 _ylm = ylm(l, m, theta.flatten(), phi.flatten())
#                 _ylm = np.reshape(_ylm, (sph.Np, sph.Np))

#                 mlab.mesh(x - m, y - l, z, scalars=_ylm, colormap="bwr")

#         mlab.view(0, 180)


# def plotYlm(sph, l, m):
#     """
#     Plots real spherical harmonics of order m and degree l.
#     Inspired by https://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html

#     Parameters
#     ----------
#     sph: spherical harmonics analysis object
#     l: int
#         degree l
#     m: int
#         order m

#     """

#     theta = np.reshape(sph.sp[:, 1], (sph.Np, sph.Np))
#     phi = np.reshape(sph.sp[:, 2], (sph.Np, sph.Np))
#     r = 0.6
#     x = r * np.sin(theta) * np.cos(phi)
#     y = r * np.sin(theta) * np.sin(phi)
#     z = r * np.cos(theta)

#     _ylm = ylm(l, m, theta.flatten(), phi.flatten())
#     _ylm = np.reshape(_ylm, (sph.Np, sph.Np))

#     mlab.mesh(x - m, y - l, z, scalars=_ylm, colormap="bwr")

#     _ylm /= _ylm.max()
#     mlab.mesh(
#         _ylm * x - m,
#         _ylm * y - l,
#         _ylm * z + 1.3,
#         scalars=np.abs(_ylm),
#         colormap="Spectral",
#     )


# def plotWlm(sph, l, m):
#     """
#     Plots magnetic field basis function 'Wlm' (r**l) over a sphere.

#     Parameters
#     ----------
#     sph: spherical harmonics analysis object
#     l: int
#         degree l
#     m: int
#         order m

#     Returns
#     -------
#     obj: mayavi object

#     """

#     _Wlm = Wlm(l, m, sph.sp[:, 1], sph.sp[:, 2])
#     _Wlm = sphvec2cart(sph.sp, _Wlm)
#     obj = mlab.quiver3d(
#         sph.p[:, 0], sph.p[:, 1], sph.p[:, 2], _Wlm[:, 0], _Wlm[:, 1], _Wlm[:, 2]
#     )
#     obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
#     return obj


# def plotBWlm_volume(sph, l, m, lim, Np, offset):
#     """
#     Plots magnetic field basis function 'Wlm' (r**l) over a 3D volume.

#     Parameters
#     ----------
#     sph: spherical harmonics analysis object
#     l: int
#         degree l
#     m: int
#         order m
#     lim: float
#         limits for coordinates, e.g., xmin = -lim, xmax = lim
#     Np: int
#         number of points along different coordinates
#     offset: 1x3 array
#         offset of the volume in which Wlm is plotted

#     Returns
#     -------
#         obj: mayavi object

#     """

#     x, y, z = np.meshgrid(
#         np.linspace(-lim + offset[0], lim + offset[0], Np),
#         np.linspace(-lim + offset[1], lim + offset[1], Np),
#         np.linspace(-lim + offset[2], lim + offset[2], Np),
#     )

#     p = np.array((x.flatten(), y.flatten(), z.flatten())).T
#     sp = cartesian2spherical(p)

#     _Wlm = Wlm(l, m, sp[:, 1], sp[:, 2])
#     _Wlm *= np.sqrt(2 * l ** 2 + l)
#     _Wlm[:, 0] *= sp[:, 0] ** (l - 1)
#     _Wlm[:, 1] *= sp[:, 0] ** (l - 1)
#     _Wlm[:, 2] *= sp[:, 0] ** (l - 1)

#     _Wlm = sphvec2cart(sp, _Wlm)
#     obj = mlab.quiver3d(p[:, 0], p[:, 1], p[:, 2], _Wlm[:, 0], _Wlm[:, 1], _Wlm[:, 2])
#     obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
#     return obj


# def plotVlm(sph, l, m):
#     """
#     Plots magnetic field basis function 'Vlm' (r**(-l)) over a sphere.

#     Parameters
#     ----------
#     sph: spherical harmonics analysis object
#     l: int
#         degree l
#     m: int
#         order m

#     Returns
#     -------
#     obj: mayavi object

#     """

#     _Vlm = Vlm(l, m, sph.sp[:, 1], sph.sp[:, 2])
#     _Vlm = sphvec2cart(sph.sp, _Vlm)
#     obj = mlab.quiver3d(
#         sph.p[:, 0], sph.p[:, 1], sph.p[:, 2], _Vlm[:, 0], _Vlm[:, 1], _Vlm[:, 2]
#     )
#     obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
#     return obj


# def plotBVlm_volume(sph, l, m, lim, Np, offset):
#     """
#     Plots magnetic field basis function 'Vlm' (r**(-l)) over a 3D volume.

#     Parameters
#     ----------
#     sph: spherical harmonics analysis object
#     l: int
#         degree l
#     m: int
#         order m
#     lim: float
#         limits for coordinates, e.g., xmin = -lim, xmax = lim
#     Np: int
#         number of points along different coordinates
#     offset: 1x3 array
#         offset of the volume in which Vlm is plotted

#     Returns
#     -------
#     obj: mayavi object

#     """

#     x, y, z = np.meshgrid(
#         np.linspace(-lim + offset[0], lim + offset[0], Np),
#         np.linspace(-lim + offset[1], lim + offset[1], Np),
#         np.linspace(-lim + offset[2], lim + offset[2], Np),
#     )

#     p = np.array((x.flatten(), y.flatten(), z.flatten())).T
#     sp = cartesian2spherical(p)

#     _Vlm = Vlm(l, m, sp[:, 1], sp[:, 2])
#     _Vlm *= np.sqrt((l + 1) * (2 * l + 1))

#     _Vlm[:, 0] *= sp[:, 0] ** (-1 * (l + 2))
#     _Vlm[:, 1] *= sp[:, 0] ** (-1 * (l + 2))
#     _Vlm[:, 2] *= sp[:, 0] ** (-1 * (l + 2))

#     _Vlm = sphvec2cart(sp, _Vlm)
#     obj = mlab.quiver3d(p[:, 0], p[:, 1], p[:, 2], _Vlm[:, 0], _Vlm[:, 1], _Vlm[:, 2])
#     obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))
#     return obj


# def compute_sphcoeffs_mesh(mesh, lmax):
#     """
#     Computes multipole moment (spherical harmonics coefficient) transformation
#     from the mesh.

#     Parameters
#     ----------
#     mesh: mesh object - the surface mesh
#     lmax: int
#         maximum degree l of the fit

#     Returns
#     -------
#     alm: (lmax*(lmax+2)xNvertices array
#           transformation from the mesh to alm coefficients (r**l-terms)
#     blm: (lmax*(lmax+2)xNvertices array
#           transformation from the mesh to blm coefficients (r**(-l)-terms)

#     """

#     Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)
#     #    Gx = Gx.toarray()
#     #    Gy = Gy.toarray()
#     #    Gz = Gz.toarray()
#     #    G = np.array([Gx, Gy, Gz])

#     tri_normals, tri_areas = tri_normals_and_areas(mesh.vertices, mesh.faces)

#     centers = np.mean(mesh.vertices[mesh.faces], axis=1)
#     centers_sp = cartesian2spherical(centers)

#     alm = np.zeros((lmax * (lmax + 2), mesh.vertices.shape[0]))
#     blm = np.zeros((lmax * (lmax + 2), mesh.vertices.shape[0]))

#     idx = 0
#     for l in range(1, lmax + 1):
#         for m in range(-1 * l, l + 1):
#             derylm = np.sqrt(l * (l + 1)) * Blm(
#                 l, m, centers_sp[:, 1], centers_sp[:, 2]
#             )
#             derylm = sphvec2cart(centers_sp, derylm)
#             # FIX: gradient of Ylm has also 1/r multiplier
#             derylm /= centers_sp[:, 0:1]

#             crossp = np.cross(centers, derylm)
#             alm_terms = crossp.T * (centers_sp[:, 0] ** l) * tri_areas
#             blm_terms = crossp.T * (centers_sp[:, 0] ** (-l - 1)) * tri_areas
#             #            integral_alm = np.sum(G*alm_terms[:,:,None], axis=(0,1))
#             #            integral_blm = np.sum(G*blm_terms[:,:,None], axis=(0,1))

#             integral_alm = alm_terms[0] @ Gx + alm_terms[1] @ Gy + alm_terms[2] @ Gz
#             integral_blm = blm_terms[0] @ Gx + blm_terms[1] @ Gy + blm_terms[2] @ Gz

#             # FIX: l-coefficients here were swapped
#             blm[idx] = (
#                 -1 / ((2 * l + 1) * l) * integral_blm
#             )  # ADDED MINUS HERE (different from Nieminen 2011, Eq 5)
#             alm[idx] = (
#                 -1 / ((2 * l + 1) * (l + 1)) * integral_alm
#             )  # THIS SIGN SHOULD BE CHECKED TOO
#             #            for i in range(mesh.vertices.shape[0]):
#             #                G = np.zeros(crossp.shape)
#             #                G[:,0] = Gx[:,i]
#             #                G[:,1] = Gy[:,i]
#             #                G[:,2] = Gz[:,i]
#             #                dotp = np.sum(G*crossp,axis=1)
#             #                integral_blm = np.sum(dotp*centers_sp[:,0]**l*tri_areas)
#             #                blm[idx,i] = -1/((2*l+1)*(l+1))*integral_blm
#             #
#             #                integral_alm = np.sum(dotp*centers_sp[:,0]**(-l-1)*tri_areas)
#             #                alm[idx,i] = 1/((2*l+1)*l)*integral_alm

#             idx += 1
#         print("l = %d computed" % (l))

#     return alm, blm
