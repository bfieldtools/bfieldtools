# -*- coding: utf-8 -*-
"""
Test sensors.py
"""

from bfieldtools.sensors import *
import numpy as np

def test_loop_magnetometer_reciprocity():
    """
    Test MagnetometerLoop 

    """
    
    m = MagnetometerLoop((0.01, 0.01), quad_scheme=("dunavant", 3))

    m1 = m.new_from_transform(np.eye(4))
    t2 = np.eye(4)
    t2[:3, 3] = np.array([0.05, 0, 0])
    m2 = m.new_from_transform(t2)

    y1 = m2.measure_bfield(m1.bfield_self)
    y2 = m2.measure_afield(m1.afield_self, ("gauss_legendre", 2))
    
    # Both measurements should yield roughly the same flux
    assert abs((y1 - y2) / y1) < 1e-2
    
    
def test_loop_gradiometer_reciprocity():
    """
    Test GradiometerLoop 

    """
    
    m = GradiometerLoop((0.01, 0.01), 0.01, quad_scheme=("dunavant", 3))

    m1 = m.new_from_transform(np.eye(4))
    t2 = np.eye(4)
    t2[:3, 3] = np.array([0.05, 0, 0])
    m2 = m.new_from_transform(t2)

    y1 = m2.measure_bfield(m1.bfield_self)
    y2 = m2.measure_afield(m1.afield_self, ("gauss_legendre", 2))

    # Both measurements should yield roughly the same flux
    assert abs((y1 - y2) / y1) < 1e-2


def test_sensor_arrays():
    points = np.array([[0, 0, 0,], [0.1, 0, 0],])

    for func in (create_mag102, create_grad204, create_arr306):
        arr = func()
        b = arr.bfields_self(points)
        a = arr.afields_self(points)
        p = arr.integration_points

        m = MagnetometerLoop((0.01, 0.01))
        m.finalize()
        f = arr.fluxes(m.bfield_self, field_type="B")
        f = arr.fluxes(m.afield_self, field_type="A")
        
        
def test_plots():
    m = MagnetometerLoop((0.1, 0.1))
    m.finalize()
    m.plot_quad_points()
    m.plot_loop()
    
    
    g = GradiometerLoop((0.05, 0.05), 0.1)
    g.finalize()
    g.plot_quad_points()
    g.plot_loop()
    


def test_measure_afield_constant_single_segment():
    sensor = MagnetometerLoop((1.0, 1.0))
    sensor.finalize()

    def afield(points):
        return np.tile([1.0, 2.0, 3.0], (len(points), 1))

    flux = sensor.measure_afield(
        afield,
        method=["gauss_legendre", 1],
    )

    np.testing.assert_allclose(flux, 0.0, atol=1e-14)
    
    
def test_stokes_theorem_unit_bfield():
    sensor = MagnetometerLoop((2.0, 3.0))
    sensor.finalize()

    def afield(points):
        out = np.empty_like(points)
        out[:, 0] = -0.5 * points[:, 1]
        out[:, 1] =  0.5 * points[:, 0]
        out[:, 2] = 0.0
        return out

    flux = sensor.measure_afield(
        afield,
        method=["gauss_legendre", 1],
    )

    np.testing.assert_allclose(
        flux,
        sensor.area,
        rtol=1e-14,
        atol=1e-14,
    )
    
    
    
def test_measure_afield_stokes():
    sensor = MagnetometerLoop((2.0, 3.0))
    sensor.finalize()

    def afield(points):
        out = np.empty_like(points)
        out[:, 0] = -0.5 * points[:, 1]
        out[:, 1] = 0.5 * points[:, 0]
        out[:, 2] = 0.0
        return out

    flux = sensor.measure_afield(
        afield,
        method=["gauss_legendre", 1],
    )

    np.testing.assert_allclose(
        flux,
        sensor.area,
        rtol=1e-12,
        atol=1e-12,
    )


def test_measure_afield_quadrature_convergence():
    m = MagnetometerLoop((0.02, 0.02))
    m1 = m.new_from_transform(np.eye(4))

    t = np.eye(4)
    t[:3, 3] = [0.05, 0, 0]
    m2 = m.new_from_transform(t)

    gl1 = m2.measure_afield(
        m1.afield_self,
        method=["gauss_legendre", 1],
    )

    gl2 = m2.measure_afield(
        m1.afield_self,
        method=["gauss_legendre", 2],
    )

    np.testing.assert_allclose(gl1, gl2, rtol=1e-3)
    
def test_magnetometer_quad_centroid():
    sensor = MagnetometerLoop(
        (2.0, 3.0),
        quad_scheme=("dunavant", 1),
    )
    sensor.finalize()

    centroid = np.average(
        sensor.area_points,
        axis=0,
        weights=sensor.quad_weights,
    )

    np.testing.assert_allclose(
        centroid,
        [0.0, 0.0, 0.0],
        atol=1e-12,
    )
    
    
def test_gradiometer_quad_centroid_total():
    sensor = GradiometerLoop(
        dimensions=(0.007, 0.021),
        baseline=0.014,
        quad_scheme=("dunavant", 3),
    )
    sensor.finalize()

    centroid = np.average(
        sensor.area_points,
        axis=0,
        weights=np.abs(sensor.quad_weights),
    )

    np.testing.assert_allclose(
        centroid,
        [0.0, 0.0, 0.0],
        atol=1e-12,
    )



