"""Loop integration, independent transient solutions, and bundled geometry."""

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.constants import mu_0
from scipy.linalg import eigh

from bfieldtools.volume import (
    TetraTopology,
    TreeCotreeGauge,
    build_edge_current_basis,
    load_example_tetra_mesh,
    loop_coupling,
    modal_response,
    reconstruct_currents,
    tetrahedron_quadrature,
)
from .test_volume_topology import cube_mesh


def basis_for(mesh):
    topology = TetraTopology.from_mesh(mesh)
    return build_edge_current_basis(
        mesh, topology, TreeCotreeGauge.from_topology(topology)
    )


def loop():
    return np.array(
        [
            [-1.0, -1.0, 2.0],
            [2.0, -1.0, 2.0],
            [2.0, 2.0, 2.0],
            [-1.0, 2.0, 2.0],
            [-1.0, -1.0, 2.0],
        ]
    )


def test_stored_mesh():
    mesh = load_example_tetra_mesh()
    assert len(mesh.tetrahedra) == 808
    assert len(mesh.vertices) == 226
    assert np.max(np.linalg.norm(mesh.vertices, axis=1)) == pytest.approx(1)
    assert mesh.cell_quality.min() > 0.25
    assert basis_for(mesh).n_dofs == 649
    assert not mesh.vertices.flags.writeable
    with pytest.raises(ValueError):
        load_example_tetra_mesh("missing")


def test_loop_coupling_against_independent_double_quadrature():
    basis = basis_for(cube_mesh(n=1))
    mesh = basis.mesh
    bary, weights = tetrahedron_quadrature(6)
    points = np.einsum("qv,cvi->cqi", bary, mesh.vertices[mesh.tetrahedra])
    u, w = leggauss(32)
    potential = np.zeros_like(points)
    wire = loop()
    for a, b in zip(wire[:-1], wire[1:]):
        samples = a + (u[:, None] + 1) / 2 * (b - a)
        inverse = 1 / np.linalg.norm(points[:, :, None, :] - samples, axis=-1)
        potential += (inverse @ (w / 2))[:, :, None] * (b - a)
    potential *= mu_0 / (4 * np.pi)
    current = reconstruct_currents(basis, np.ones(basis.n_dofs))
    reference = np.einsum("cqi,ci,q,c->", potential, current, weights, mesh.volumes)
    actual = loop_coupling(basis, wire, quadrature_order=6, cell_chunk_size=2)
    assert actual.sum() == pytest.approx(reference, rel=1e-9, abs=1e-20)
    np.testing.assert_allclose(
        loop_coupling(basis, wire[::-1], quadrature_order=6),
        -actual,
        rtol=1e-10,
        atol=1e-20,
    )


@pytest.mark.parametrize(
    "wire",
    [
        [[-1, 0.5, 0.5], [2, 0.5, 0.5], [2, 2, 2], [-1, 0.5, 0.5]],
        [[0, 0, 0], [-1, 0, 0], [-1, -1, 0], [0, 0, 0]],
        [[0, 0, 2], [1, 0, 2], [0, 1, 2]],
        [[0, 0, 2], [1, 0, 2], [1, 0, 2], [0, 0, 2]],
    ],
)
def test_invalid_or_intersecting_loops(wire):
    with pytest.raises(ValueError):
        loop_coupling(basis_for(cube_mesh(n=1)), wire)


def test_modal_response_constant_ramp_and_free_decay():
    rates = np.array([1.0, 2.0, 10.0])
    coupling = np.array([0.5, -0.3, 0.2])
    times = np.linspace(0, 2, 81)
    initial = np.array([0.2, -0.1, 0.7])
    result = modal_response(
        rates, coupling, lambda t: 3.0, times, initial=initial, rtol=1e-11
    )
    expected = initial[:, None] * np.exp(-rates[:, None] * times) - 3 * coupling[
        :, None
    ] / rates[:, None] * (1 - np.exp(-rates[:, None] * times))
    np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-10)
    free = modal_response(rates, coupling, lambda t: 0.0, times, initial=initial)
    np.testing.assert_allclose(
        free, initial[:, None] * np.exp(-rates[:, None] * times), rtol=2e-5, atol=1e-11
    )


def test_modal_response_sign_metric_and_energy():
    resistance = np.array([[3.0, 0.4], [0.4, 2.0]])
    magnetic = np.array([[2.0, 0.2], [0.2, 1.0]])
    coupling = np.array([0.3, -0.1])
    rates, vectors = eigh(resistance, magnetic)
    times = np.linspace(0, 1, 50)
    a = modal_response(rates, vectors.T @ coupling, lambda t: 1.0, times)
    x = vectors @ a
    da = -rates[:, None] * a - (vectors.T @ coupling)[:, None]
    dx = vectors @ da
    np.testing.assert_allclose(
        magnetic @ dx + resistance @ x,
        -np.broadcast_to(coupling[:, None], x.shape),
        atol=1e-14,
    )
    # d(1/2 x.T M x)/dt = -x.T R x - x.T m dI/dt.
    np.testing.assert_allclose(
        np.sum(x * (magnetic @ dx), axis=0),
        -np.sum(x * (resistance @ x), axis=0) - coupling @ x,
        atol=1e-14,
    )


@pytest.mark.parametrize(
    "rates,coupling,times",
    [
        ([-1], [1], [0, 1]),
        ([1], [1, 2], [0, 1]),
        ([1], [1], [1, 0]),
        ([1], [1], [0, np.nan]),
        ([1j], [1], [0, 1]),
    ],
)
def test_invalid_modal_inputs(rates, coupling, times):
    with pytest.raises(ValueError):
        modal_response(rates, coupling, lambda t: 0.0, times)
