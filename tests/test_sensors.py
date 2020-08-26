# -*- coding: utf-8 -*-
"""
Test sensors.py
"""

from bfieldtools.sensors import *


def test_loop_sensors():
    """
    Test MagnetometerLoop and GradiometerLoop classes

    """
    for cl in MagnetometerLoop, GradiometerLoop:
        m = cl((1, 1))

        m1 = m.new_from_transform(np.eye(4))
        t2 = np.eye(4)
        t2[:3, 3] = np.array([3, 0, 0])
        m2 = m.new_from_transform(t2)

        points = np.array([[1, 0, 1,], [2, 0, 0]])

        b = m1.bfield_self(points)
        a = m1.afield_self(points)

        y1 = m2.measure_bfield(m1.bfield_self)
        y2 = m2.measure_afield(m1.afield_self)

        # Both measurements should yield roughly the same flux
        assert abs((y1 - y2) / y1) < 1e-2

        A = m2.area
        p = m2.integration_points

        # TODO test functionality with analytic stuff


def test_sensor_arrays():
    points = np.array([[0, 0, 0,], [0.1, 0, 0]])

    for func in (create_mag102, create_grad204, create_arr306):
        arr = func()
        b = arr.bfields_self(points)
        a = arr.afields_self(points)
        p = arr.integration_points
