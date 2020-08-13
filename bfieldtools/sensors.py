# -*- coding: utf-8 -*-
"""
Classes for magnetic sensors and sensor arrays
"""

import quadpy
import numpy as np
from abc import ABC, abstractmethod
from .line_conductor import LineConductor


def get_default_scheme(geometry, Nq):
    if isinstance(Nq, tuple) or isinstance(Nq, list):
        Nq = "_".join(map(str, Nq))
    if geometry == "c2":
        if Nq is None:
            Nq = "01"
        return quadpy.c2.__all__[f"dunavant_{Nq}"]()
    if geometry == "c3":
        return quadpy.c3.__all__[f"hammer_stroud_{Nq}"]()


# How to define an abstract base class
class BaseSensor(ABC):
    def __init__(self):
        super().__init__()

    def new_from_transform(self, transform):
        obj = self.__class__(**self.__dict__)
        obj.finalize(transform)

        return obj

    @abstractmethod
    def finalize(self, transform):
        t = transform
        if not isinstance(t, np.ndarray):
            raise TypeError("transform must be ndarray")
        if not t.shape == (4, 4):
            raise ValueError("transform shape must be (4,4)")

    @abstractmethod
    def update_position(self, transform):
        t = transform
        if not isinstance(t, np.ndarray):
            raise TypeError("transform must be ndarray")
        if not t.shape == (4, 4):
            raise ValueError("transform shape must be (4,4)")

    @abstractmethod
    def quad_points(self,):
        pass

    @abstractmethod
    def bfield(self, points):
        pass

    @abstractmethod
    def afield(self, points):
        pass


class MagnetometerLoop(BaseSensor):
    def __init__(self, dimensions, quad_scheme="dunavant_01"):
        if len(dimensions) == 2:
            self.dimensions = tuple(dimensions)
        else:
            raise ValueError("len(dimensions) must be 2")
        self.quad_scheme = quad_scheme

    def finalize(self, transform=np.eye(4)):
        super().finalize(transform)

        r0 = self.dimensions[0] / 2
        r1 = self.dimensions[1] / 2
        points = np.zeros(4, 3)
        points[:, 0] = np.array([-r0, r0, r0, -r0])
        points[:, 1] = np.array([-r1, -r1, r1, r1])
        self.points = points
        self.update_position(transform)

    def update_position(self, transform=np.eye(4)):
        super().update_position()
        points = apply_transform(transform, self.points)
        self.line_conductor = LineConductor((points))

        # TODO quad points

    def quad_points(self):
        return quadpy.c2.__all__[self.quad_scheme].points

    def integrate(self, func):
        quadpy.c2.__all__[self.quad_scheme]

    def bfield(self):
        pass

    def afield():
        pass


class GradiometerLoop(BaseSensor):
    def __init__(self, dimensions, baseline, scheme=None, Nq=None):
        pass


class MagnetometerOPM(BaseSensor):
    def __init__(self, dimensions, scheme=None, Nq=None):
        pass


def apply_transform(transform, points):
    """
    Parameters
    ----------
    transform : numpy array (4,4)
        4x4 transformation matrix for homogeneous coordinates
    points : numpy array (N, 3)
        points to be transformed

    Returns
    -------
    transformed points

    """
    return (transform[:3, :3] @ points.T + transform[:3, 3:]).T


class SensorArray:
    def __init__(self, base_sensor, transforms, Nq):
        self.transforms = transforms
        self.sensors = []

        for t in self.transforms:
            sensor = base_sensor.new_from_transform(t)
            self.sensors.append(sensor)

    def bfields(self, points):
        return np.array([s.bfield(points) for s in self.sensors])

    def afields(self, points):
        return np.array([s.afield(points) for s in self.sensors])
