"""
This module includes a convenience class for working with polyline currents
"""

__all__ = ["LineConductor"]

import numpy as np
from trimesh.path import Path3D
from trimesh.path.entities import Line

from .viz import plot_3d_current_loops
from .contour import scalar_contour, simplify_contour
from .mesh_impedance import mesh2line_mutual_inductance
from . import line_magnetics


class LineConductor(Path3D):
    """
    Class that inherits Trimesh.path.Path3D for handling discretized current loops.
    Functions inside assume that vertices are unique for each entity.

    """

    def __init__(self, loops=None, mesh=None, scalars=None, **kwargs):
        """
        Init with inheritance. First priority is given to passed loops parameter,
        if not present will compute loops from mesh and scalars.

        Parameters
        ----------
        loops: list of array-like (x, 3)
            Each list element corresponds to a Polyline (current loop), with the
            last element being connected to the first
        mesh: Trimesh mesh object
            mesh geometry
        scalars: array-like or StreamFunction
            scalar function defined in the vertices of the mesh
        kwargs: dict
            passed to scalar_contour if called. Relevant kw:s are N_contours
            and contours
            
        """
        vertices = np.zeros((0, 3))
        entities = []

        if loops is None:
            loops, vals = scalar_contour(mesh, scalars, return_values=True, **kwargs)

        for loop in loops:
            if np.all(loop[0] == loop[-1]):
                points = np.arange(0, len(loop)) + len(vertices)
            else:  # Enforce closed loops
                points = np.append(np.arange(0, len(loop)), 0) + len(vertices)
            entities.append(Line(points))

            vertices = np.append(vertices, loop, axis=0)

        Path3D.__init__(self, entities, vertices)

        if loops is None:
            self.vals = vals

    def simplify(self, min_edge=1e-3, angle_threshold=2e-2, smooth=True):
        """
        Simplifies contour paths.
        
        Parameters
        ----------
        c: array-like
            List of polygons describing closed loops.
        min_edge: float
            Minimum edge length. Edges shorter than this are merged.
        angle_threshold: float
            Minimum angle. Edges with smaller angle differences are merged.
        smooth: bool
            If True, apply smoothing to the polygon shapes.

        Returns
        -------
        simplified_linepath: LineConductor
        
        """
        simplified_loops = [
            simplify_contour(
                e.discrete(self.vertices), min_edge, angle_threshold, smooth
            )
            for e in self.entities
        ]

        return LineConductor(simplified_loops)

    def plot_loops(self, **kwargs):
        """
        Plots loops in 3D using mayavi, see viz.plot_3d_current_loops for more details

        Parameters
        ----------
        colors: str


        """

        if "tube_radius" not in kwargs:
            kwargs["tube_radius"] = 0.005 * np.linalg.norm(self.bounds)

        figure = plot_3d_current_loops(
            [e.discrete(self.vertices) for e in self.entities], **kwargs
        )
        return figure

    def magnetic_field(self, points, separate_loops=False):
        """
        Compute magnetic field in some point due to a unit current in the loops

        Parameters
        ----------
        points: array (N_p, 3)

        separate_loops: Boolean
            	If True, don't combine contributions of separate loops

        Returns
        -------

        Bfield: array (N_p, 3) or array (N_loops, N_p, 3)

        """
        Bfield = np.zeros((len(self.entities), len(points), 3))
        for ii, loop in enumerate(self.entities):
            Bfield[ii] = line_magnetics.magnetic_field(
                self.vertices[loop.points], points
            )

        if not separate_loops:
            Bfield = np.sum(Bfield, axis=0)

        return Bfield

    def vector_potential(self, points, separate_loops=False, **kwargs):
        """
        Compute magnetic vector potential in some point due to a unit current in the loops

        Parameters
        ----------
        points: array (N_p, 3)

        separate_loops: Boolean
            	If True, don't combine contributions of separate loops

        Returns
        -------

        Aield: array (N_p, 3) or array (N_loops, N_p, 3)
        """

        Afield = np.zeros((len(self.entities), len(points), 3))
        for ii, loop in enumerate(self.entities):
            Afield[ii] = line_magnetics.vector_potential(
                self.vertices[loop.points], points, **kwargs
            )

        if not separate_loops:
            Afield = np.sum(Afield, axis=0)

        return Afield

    def scalar_potential(self, points, separate_loops=False, **kwargs):
        """
        Compute magnetic scalar potential in some point due to a unit current in the loops

        Parameters
        ----------
        points: array (N_p, 3)

        separate_loops: Boolean
            	If True, don't combine contributions of separate loops

        Returns
        -------

        Ufield: array (N_p,) or array (N_loops, N_p)
        """

        Ufield = np.zeros((len(self.entities), len(points)))
        for ii, loop in enumerate(self.entities):
            Ufield[ii] = line_magnetics.scalar_potential(
                self.vertices[loop.points], points, **kwargs
            )

        if not separate_loops:
            Ufield = np.sum(Ufield, axis=0)

        return Ufield

    def line_mutual_inductance(self, path, separate_loops=False, **kwargs):
        """
        Calculate mutual inductance between self and another line current (path)

        Parameters
        ----------
        path : Path3d or LineConductor object
            The other line conductor
        separate_loops : TYPE, optional
            If True, return the inductance separately for each loop. The default is False.

        Returns
        -------
        M: array
            If separate loops, shape (N_loops, N_loops_path). Otherwise (N_loops_path,).
            where N_loops_path = len(path.entities)

        """

        M = line_magnetics.mutual_inductance(self, path, **kwargs)

        if not separate_loops:
            M = np.sum(M, axis=0)

        return M

    def mesh_mutual_inductance(self, mesh, separate_loops=False, **kwargs):
        """
        Computes the mutual inductance between the polyline loops and a mesh_conductor
        mesh.

        Parameters
        ----------
        mesh : Trimesh mesh object
            Mesh with N_verts vertices.
        separate_loops : Boolean, optional
            If True, return the inductance separately for each loop. The default is False.
        **kwargs : dict
            Passed to mesh_impedance.mesh2line_mutual_inductance

        Returns
        -------
        M: array
            If separate loops, shape (N_loops, N_verts). Otherwise (N_verts,).

        """
        M = np.zeros((len(self.entities), len(mesh.vertices)))

        for ii, loop in enumerate(self.entities):
            M[ii] = mesh2line_mutual_inductance(
                mesh, self.vertices[loop.points], **kwargs
            )

        if not separate_loops:
            M = np.sum(M, axis=0)

        return M
