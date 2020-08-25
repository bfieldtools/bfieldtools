# -*- coding: utf-8 -*-
"""
Test sensors.py
"""

from bfieldtools.sensors import *


def test_magnetometerloop():
    """
    Test MagnetometerLoop class

    """

    m = MagnetometerLoop((1, 1))
    # make new magnetometer using the following method
    m1 = m.new_from_transform(np.eye(4))
    t2 = np.eye(4)
    t2[:3, 3] = np.array([3, 0, 0])
    m2 = m.new_from_transform(t2)

    points = np.array([[1, 0, 1,], [2, 0, 0]])

    b = m1.bfield_self(points)
    a = m1.afield_self(points)

    m2.measure_bfield(m1.bfield_self)
    m2.measure_afield(m1.afield_self)

    A = m2.area
    p = m2.integration_points

    # TODO test functionality with analytic stuff


def test_sensor_array():
    arr = create_mag102()

    points = np.array([[0, 0, 0,], [0.1, 0, 0]])
    b = arr.bfields_self(points)
    a = arr.afields_self(points)
