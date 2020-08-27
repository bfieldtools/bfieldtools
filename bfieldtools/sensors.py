# -*- coding: utf-8 -*-
"""
Classes for magnetic sensors and sensor arrays
"""

import quadpy
import numpy as np
from abc import ABC, abstractmethod
from .line_conductor import LineConductor


class BaseSensor(ABC):
    def __init__(self):
        self.base_sensor = True
        super().__init__()

    def new_from_transform(self, transform):
        """
        Create a copy of the sensor and finalize it with 'transform'

        Parameters
        ----------
        transform : ndarray (4,4)
           transformation matrix


        Returns
        -------
        obj : self.__class__
            New finalized object of the same class
            and initialization parameters

        """
        if self.base_sensor:
            d = self.__dict__.copy()
            d.pop("base_sensor")
            obj = self.__class__(**d)
            obj.finalize(transform)
            return obj
        else:
            raise TypeError("Sensor already finalized, to move use update_position")

    @abstractmethod
    def finalize(self, transform):
        """
        Finalize the geometric information of the sensor
        and transform it using 'transform'

        After finalize the object is ready for field measurements
        and reciprocal field calculations

        Parameters
        ----------
        transform : ndarray (4,4)
           transformation matrix


        Returns
        -------
        None

        """
        self.base_sensor = False
        t = transform
        if not isinstance(t, np.ndarray):
            raise TypeError("transform must be ndarray")
        if not t.shape == (4, 4):
            raise ValueError("transform shape must be (4,4)")

    @abstractmethod
    def update_position(self, transform):
        """
        Update the geometric information of the sensor
        by transform

        transform acts on the initial coordinates of the sensor
        (not on the transformed coordinates)

        Parameters
        ----------
        transform : ndarray (4,4)
           transformation matrix. Transforms from sensor coordinate system
           to array coordinate system

        Returns
        -------
        None

        """
        t = transform
        if not isinstance(t, np.ndarray):
            raise TypeError("transform must be ndarray")
        if not t.shape == (4, 4):
            raise ValueError("transform shape must be (4,4)")

    @abstractmethod
    def measure_bfield(self, bfield_func):
        """
        Magnetic field measurement

        Parameters
        ----------
        bfield_func : function
            function taking (N, 3) array of points
            as parameter and outputting (N, 3) array of magnetic field

        Returns
        -------
        meas : float
            Magnetic field measurement

        """
        pass

    @abstractmethod
    def measure_afield(self, afield_func):
        """
        Magnetic field measurement using vector potential

        Parameters
        ----------
        afield_func : function
            function that takes (N, 3) array of points
            as a parameter and outputs (N, 3) array of vector potential

        Returns
        -------
        meas : float
            Magnetic field measurement using vector potential

        """
        pass

    @abstractmethod
    def bfield_self(self, points):
        """
        Reciprocal B-field for the sensor

        Parameters
        ----------
        points : ndarray (N, 3)
            evaluation points

        Returns
        -------
        B-field at the evaluation points

        """
        pass

    @abstractmethod
    def afield_self(self, points):
        """
        Reciprocal vector potential for the sensor

        Parameters
        ----------
        points : ndarray (N, 3)
            evaluation points

        Returns
        -------
        B-field at the evaluation points

        """
        pass


class MagnetometerLoop(BaseSensor):
    """
    A class for rectangular pickup-loop based magnetometers

    Initialization
    --------------

    dimensions : tuple/list of loop side lengths (float)

    quad_scheme : quadrature scheme for bfield integration
                 Default: "dunavant_01"
                 for more schemes see quadpy.c2.__all__
    """

    def __init__(self, dimensions, quad_scheme="dunavant_01"):
        if len(dimensions) == 2:
            self.dimensions = tuple(dimensions)
        else:
            raise ValueError("len(dimensions) must be 2")
        self.quad_scheme = quad_scheme
        super().__init__()

    def init_geometry(self):
        r0 = self.dimensions[0] / 2
        r1 = self.dimensions[1] / 2
        points = np.zeros((4, 3))
        points[:, 0] = np.array([-r0, r0, r0, -r0])
        points[:, 1] = np.array([-r1, -r1, r1, r1])
        self.loop_points = points

    def finalize(self, transform=np.eye(4)):
        t = transform
        super().finalize(t)

        self.init_geometry()

        self.quad_weights = self.scheme.weights
        self.area_points = np.zeros((len(self.scheme.points), 3))
        self.area_points[:, 0] = self.scheme.points[:, 0] * self.dimensions[0] / 2
        self.area_points[:, 1] = self.scheme.points[:, 1] * self.dimensions[1] / 2

        self.update_position(t)

    def update_position(self, transform=np.eye(4)):
        t = transform
        super().update_position(t)
        points = apply_transform(t, self.loop_points)
        self.line_conductor = LineConductor([points])

        self.bfield_points = apply_transform(t, self.area_points)
        self.normal = t[:3, 2] / np.linalg.norm(t[:3, 2])

    @property
    def area(self):
        return np.prod(self.dimensions)

    @property
    def integration_points(self):
        return self.bfield_points

    @property
    def scheme(self):
        return quadpy.c2.__dict__[self.quad_scheme]()

    def measure_bfield(self, bfield_func):
        bfield = bfield_func(self.bfield_points)
        integral = np.einsum("ij,i,j->", bfield, self.quad_weights, self.normal)
        return integral * self.area

    def measure_afield(self, afield_func, scheme_degree=1):
        scheme = quadpy.c1.gauss_patterson(scheme_degree)
        lc = self.line_conductor
        segments = lc.discrete[0][1:] - lc.discrete[0][:-1]
        # Points are from -1 to 1
        w = (scheme.points + 1) / 2
        quad_points = [
            p + w[:, None] * s for s, p in zip(segments, lc.discrete[0][:-1])
        ]
        # Nsegment, Nqpoints_per_segment, Nxyz
        quad_points = np.array(quad_points)
        shape = quad_points.shape
        afield = afield_func(quad_points.reshape(-1, 3)).reshape(shape)
        # scheme.weights sum to line length == 2
        weights = scheme.weights / 2

        # Sum integral using einsum
        return np.einsum("ijk,ik,j", afield, segments, weights)

    def bfield_self(self, points):
        return self.line_conductor.magnetic_field(points)

    def afield_self(self, points):
        return self.line_conductor.vector_potential(points)

    def plot_loop(self, color=(0, 0, 1)):
        self.line_conductor.plot_loops(figure=False, tube_radius=None, colors=color)

    def plot(self):
        self.plot_loop()

    def plot_quad_points(self):
        from mayavi import mlab

        mlab.points3d(
            *self.bfield_points.T, color=(1, 0, 0), scale_factor=self.dimensions[0] / 5
        )


class GradiometerLoop(MagnetometerLoop):
    """
    A class for rectangular pickup-loop based gradiometers

    Initialization
    --------------

    dimensions : tuple/list of loop side lengths (float)

    baseline : float, distance of the loop centers in x direction

    quad_scheme : quadrature scheme for bfield integration
                 Default: "dunavant_01"
                 for more schemes see quadpy.c2.__all__
    """

    def __init__(
        self, dimensions=(0.007, 0.021), baseline=0.014, quad_scheme="dunavant_01"
    ):
        self.baseline = baseline
        super().__init__(dimensions, quad_scheme)

    def init_geometry(self):
        """
        Loop points for gradiometer measuring x-derivative

        """
        r0 = self.dimensions[0] / 2
        r1 = self.dimensions[1] / 2
        points = np.zeros((6, 3))
        points[:, 0] = np.array([-r0, -r0, r0, r0, -r0, -r0])
        points[:, 1] = np.array([0, -r1, -r1, r1, r1, 0])

        points2 = points.copy()
        points2[:, 0] *= -1

        points[:, 0] += self.baseline / 2
        points2[:, 0] -= self.baseline / 2

        self.loop_points = np.vstack((points, points2))

    def finalize(self, transform=np.eye(4)):
        t = transform
        BaseSensor.finalize(self, t)

        self.init_geometry()

        w = self.scheme.weights
        self.quad_weights = np.hstack((w, -w))
        area_points = np.zeros((len(self.scheme.points), 3))
        area_points[:, 0] = self.scheme.points[:, 0] * self.dimensions[0] / 2
        area_points[:, 1] = self.scheme.points[:, 1] * self.dimensions[1] / 2

        area_points2 = area_points.copy()
        area_points2[:, 0] *= -1
        area_points[:, 0] += self.baseline / 2
        area_points2[:, 0] -= self.baseline / 2

        self.area_points = np.vstack((area_points, area_points2))

        self.update_position(t)

    def plot_quad_points(self):
        from mayavi import mlab

        mask = self.quad_weights > 0
        mlab.points3d(
            *self.bfield_points[mask].T,
            color=(1, 0, 0),
            scale_factor=self.dimensions[0] / 5,
        )
        mlab.points3d(
            *self.bfield_points[~mask].T,
            color=(0, 0, 1),
            scale_factor=self.dimensions[0] / 5,
        )


class MagnetometerOPM(BaseSensor):
    """OPM sensor modelled as a point"""

    def __init__(self, axes=(0,)):
        self.axes = axes

    def finalize(transform):
        pass


class ThreeAxis(BaseSensor):
    def __init__(self, mag_points, mag_dirs=np.eye(3)):
        self.points = mag_points
        self.dirs = mag_dirs

    def finalize(self, transform):
        t = transform
        super().finalize(t)

    def update_position(self, transform):
        t = transform
        super().update_position(t)


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
    """A class for storing an array of sensors

    Initialization
    --------------
    base_sensor : a sensor object (child class of BaseSensor)
                  that is copied for all the sensor positions
    transforms : 4x4 transformation matrices for all the sensors
                these transforms map the sensor geometry from the
                sensor coordinate system to the array coordinate system
    names : list of str. Names of the sensor channels

    """

    def __init__(self, base_sensor, transforms, names):
        self.transforms = transforms
        self.names = names
        self.sensors = {}

        for k, t in zip(self.names, self.transforms):
            sensor = base_sensor.new_from_transform(t)
            self.sensors[k] = sensor

        self._integration_points = None

    def __repr__(self):
        sensor_type = list(self.sensors.values())[0].__class__.__name__
        count = len(self.names)
        return f"SensorArray: {count} sensors of type {sensor_type}"

    @property
    def integration_points(self):
        if self._integration_points is None:
            p = np.array([s.integration_points for s in self.sensors.values()])
            self._integration_points = p
        return self._integration_points

    def fluxes(self, field_func, field_type="B"):
        """
        Measure the flux either based on the magnetic field (V)
        or the vector potential (A)

        Parameters
        ----------
        field_func : function
            field function that takes points as (N, 3) ndarray as input
        field_type : str, optional
            'A' vector potential, 'B' magnetic field. The default is 'B'.
        Returns
        -------
        ndarray
            ndarray of fluxes in each sensor

        """
        if field_type == "A":
            return np.array(
                [s.measure_afield(field_func) for s in self.sensors.values()]
            )
        elif field_type == "B":
            return np.array(
                [s.measure_bfield(field_func) for s in self.sensors.values()]
            )
        else:
            raise ValueError(
                "Field type must be either A (vector potential) or B (magnetic field)"
            )

    def bfields_self(self, points):
        return np.array([s.bfield_self(points) for s in self.sensors.values()])

    def afields_self(self, points):
        return np.array([s.afield_self(points) for s in self.sensors.values()])

    def plot(self):
        for k in self.sensors.keys():
            self.sensors[k].plot()


def call_subarrays(f):
    name = f.__name__

    def call_them(*args, **kwargs):
        obj = args[0]
        args = args[1:]
        outlist = [arr.__getattribute__(name) for arr in obj.arrays]
        # call with arguments if the attribute is a function
        if hasattr(outlist[0], "__call__"):
            outlist = [func(*args, **kwargs) for func in outlist]
        if obj.concatenate and isinstance(outlist[0], np.ndarray):
            try:
                return np.concatenate(outlist, axis=0)
            except ValueError:
                print("Cannot concatenate, returning list of objects")
                return outlist
        else:
            return outlist

    return call_them


class MixedArray:
    """Array of mixed types of sensors

    Initialization
    --------------
        arrays : tuple/list of SensorArray-objects (subarrays)

        concatenate : boolean, concatenate output arrays (see below)

    The methods of this class call the subarrays
    in stored in 'arrays' using the decorator 'call_subarrays'

    If the method output is a numpy ndarray, the outputs
    can be automatically concetenated by having concatenate=True
    """

    def __init__(self, arrays, concatenate=True):
        self.arrays = arrays
        self.concatenate = concatenate

    @property
    @call_subarrays
    def integration_points(self):
        pass

    @call_subarrays
    def fluxes(self, field_func, field_type="B"):
        pass

    @call_subarrays
    def bfields_self(self, points):
        pass

    @call_subarrays
    def afields_self(self, points):
        pass

    @call_subarrays
    def plot(self):
        pass


def create_mag102():
    import pkg_resources

    fname = (
        pkg_resources.resource_filename("bfieldtools", "sensor_data/mag102_trans.npz"),
    )
    arr_info = np.load(fname[0])

    names = list(arr_info.keys())
    mats = list(arr_info.values())

    bs = MagnetometerLoop((0.021, 0.021))

    return SensorArray(bs, mats, names)


def create_grad204():
    import pkg_resources

    fname = (
        pkg_resources.resource_filename("bfieldtools", "sensor_data/grad204_trans.npz"),
    )
    arr_info = np.load(fname[0])

    names = list(arr_info.keys())
    mats = list(arr_info.values())

    bs = GradiometerLoop((0.007, 0.021), 0.014)

    return SensorArray(bs, mats, names)


def create_arr306():
    mags = create_mag102()
    grads = create_grad204()

    return MixedArray((mags, grads))
