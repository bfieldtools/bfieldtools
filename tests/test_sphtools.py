from bfieldtools import sphtools

import pytest


import numpy as np

from numpy.testing import assert_allclose


def test_coord_changes():

    point = np.array([[0.4, 0.4, 1]])

    sphpoint = sphtools.cartesian2spherical(point)

    point2 = sphtools.spherical2cartesian(sphpoint)

    assert_allclose(point, point2)


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

    p = np.array([[1, 0, 0], [-1, 2, 0]])

    # make a homogeneous field
    acoeffs = np.array([0, 0, 0, 0, 0, 0, 0, 0])
    bcoeffs = np.array([1, 0, 0, 0, 0, 0, 0, 0])
    lmax = 2

    B = sphtools.field(p, acoeffs, bcoeffs, lmax)

    assert_allclose(B[0], B[1], atol=1e-16)

    basis = sphtools.basis_fields(p, lmax)

    B_b = basis[0] @ acoeffs + basis[1] @ bcoeffs

    assert_allclose(B_b[0], B_b[1], atol=1e-16)

    assert_allclose(B, B_b, atol=1e-16)

    p = np.array([[0, 1, 0], [100, 2, 500]])

    # make a homogeneous field
    acoeffs = np.array([0, 0, 0, 0, 0, 0, 0, 0])
    bcoeffs = np.array([0, 1, 0, 0, 0, 0, 0, 0])
    lmax = 2

    B = sphtools.field(p, acoeffs, bcoeffs, lmax)

    assert_allclose(B[0], B[1], atol=1e-16)

    basis = sphtools.basis_fields(p, lmax)

    B_b = basis[0] @ acoeffs + basis[1] @ bcoeffs

    assert_allclose(B_b[0], B_b[1], atol=1e-16)

    assert_allclose(B, B_b, atol=1e-16)


def test_innerproduct():

    sph = sphtools.SphBasis(40)

    atol = 1e-12

    Ylm1 = sphtools.ylm(2, 1, sph.sqp[:, 1], sph.sqp[:, 2])
    Ylm2 = sphtools.ylm(2, 2, sph.sqp[:, 1], sph.sqp[:, 2])

    assert_allclose(sph.innerproduct(Ylm1[:, None], Ylm1[:, None]), 1, atol=atol)
    assert_allclose(sph.innerproduct(Ylm2[:, None], Ylm2[:, None]), 1, atol=atol)
    assert_allclose(sph.innerproduct(Ylm1[:, None], Ylm2[:, None]), 0, atol=atol)

    Vlm1 = sphtools.Vlm(2, 1, sph.sqp[:, 1], sph.sqp[:, 2])
    Vlm2 = sphtools.Vlm(5, 2, sph.sqp[:, 1], sph.sqp[:, 2])

    assert_allclose(sph.innerproduct(Vlm1, Vlm1), 1, atol=atol)
    assert_allclose(sph.innerproduct(Vlm2, Vlm2), 1, atol=atol)
    assert_allclose(sph.innerproduct(Vlm1, Vlm2), 0, atol=atol)

    Wlm1 = sphtools.Wlm(3, 1, sph.sqp[:, 1], sph.sqp[:, 2])
    Wlm2 = sphtools.Wlm(5, 2, sph.sqp[:, 1], sph.sqp[:, 2])

    assert_allclose(sph.innerproduct(Wlm1, Wlm1), 1, atol=atol)
    assert_allclose(sph.innerproduct(Wlm2, Wlm2), 1, atol=atol)
    assert_allclose(sph.innerproduct(Wlm1, Wlm2), 0, atol=atol)

    Xlm1 = sphtools.Xlm(3, 1, sph.sqp[:, 1], sph.sqp[:, 2])
    Xlm2 = sphtools.Xlm(4, 2, sph.sqp[:, 1], sph.sqp[:, 2])

    assert_allclose(sph.innerproduct(Xlm1, Xlm1), 1, atol=atol)
    assert_allclose(sph.innerproduct(Xlm2, Xlm2), 1, atol=atol)
    assert_allclose(sph.innerproduct(Xlm1, Xlm2), 0, atol=atol)

    assert_allclose(sph.innerproduct(Xlm1, Vlm1), 0, atol=atol)
    assert_allclose(sph.innerproduct(Xlm1, Vlm2), 0, atol=atol)
    assert_allclose(sph.innerproduct(Xlm2, Vlm1), 0, atol=atol)
    assert_allclose(sph.innerproduct(Xlm2, Vlm2), 0, atol=atol)

    assert_allclose(sph.innerproduct(Wlm1, Vlm1), 0, atol=atol)
    assert_allclose(sph.innerproduct(Wlm1, Vlm2), 0, atol=atol)
    assert_allclose(sph.innerproduct(Wlm2, Vlm1), 0, atol=atol)
    assert_allclose(sph.innerproduct(Wlm2, Vlm2), 0, atol=atol)


def test_mesh_coupling():
    """
        Test compute_sphcoeffs_mesh with sphere
    """

    from bfieldtools.sphtools import compute_sphcoeffs_mesh
    from bfieldtools.utils import load_example_mesh
    from bfieldtools.sphtools import basis_potentials, basis_fields
    from bfieldtools.sphtools import ylm, cartesian2spherical
    from bfieldtools.mesh_calculus import mass_matrix

    mesh = load_example_mesh("unit_sphere")
    mesh.vertices *= 1 / np.sqrt(mesh.area / (4 * np.pi))  # Scale to unit radius
    R = 2
    mesh.vertices *= R  # Scale to R

    c = compute_sphcoeffs_mesh(mesh, 3)

    # Test potential
    sp = cartesian2spherical(mesh.vertices)
    M = mass_matrix(mesh, lumped=True)
    u1, u2 = basis_potentials(mesh.vertices, 3)
    diff1 = []
    diff2 = []

    for ll in range(1, 4):
        for m in range(-ll, ll + 1):
            s = ylm(ll, m, sp[:, 1], sp[:, 2])
            p = u1 @ c[0] @ s
            # p should be p= ll/(2*ll+1)s, test this
            coeff = s @ M @ p / (s @ M @ s)
            diff1.append(coeff - ll / (2 * ll + 1))
            p = u2 @ c[1] @ s
            # p should be p= -(ll+1)/(2*ll+1)s, test this
            coeff = s @ M @ p / (s @ M @ s)
            diff2.append(coeff + (ll + 1) / (2 * ll + 1))

    # The integration accuracy is quite low so set the tolerance high
    assert np.allclose(diff1, 0, atol=1e-2)
    assert np.allclose(diff2, 0, atol=1e-2)

    # Test field
    b1, b2 = basis_fields(mesh.vertices, 3)
    b1 = np.einsum("ijk,ij->ik", b1, mesh.vertex_normals)
    b2 = np.einsum("ijk,ij->ik", b2, mesh.vertex_normals)
    diff1 = []
    diff2 = []
    mu0 = 4 * np.pi * 1e-7
    for ll in range(1, 4):
        for m in range(-ll, ll + 1):
            s = ylm(ll, m, sp[:, 1], sp[:, 2])
            p = b1 @ c[0] @ s
            # p should be p= mu0*(ll+1)*ll/(2*ll+1)/R s, test this
            coeff = s @ M @ p / (s @ M @ s)
            diff1.append(coeff / mu0 - (ll + 1) * ll / (2 * ll + 1) / R)
            p = b2 @ c[1] @ s
            # p should be p= mu0*(ll+1)*ll/(2*ll+1)/R s, test this
            coeff = s @ M @ p / (s @ M @ s)
            diff2.append(coeff / mu0 - (ll + 1) * ll / (2 * ll + 1) / R)

    assert np.allclose(diff1, 0, atol=1e-2)
    assert np.allclose(diff2, 0, atol=1e-2)
