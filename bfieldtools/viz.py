"""
Visualization functions tailored for bfieldtools.
Mainly wrappers and convenience helpers around pyvista and matplotlib functions
"""

__all__ = [
    "plot_3d_current_loops",
    "plot_cross_section",
    "plot_data_on_faces",
    "plot_data_on_vertices",
    "plot_mesh",
]

import matplotlib.pyplot as plt
from matplotlib import colors as c

import numpy as np


def plot_mesh(
    mesh, cull_front=False, cull_back=False, figure=None, figsize=(800, 800), **kwargs
):
    """
    Plot the mesh surface in pyvista

    Parameters
    ----------
    mesh: Trimesh mesh object
        mesh to be plotted
    cull_front: Boolean (False)
        If True, cull front of mesh
    cull_back: Boolean (False)
        If True, cull back of mesh
    figure: pyvista figure (Plotter) or None
        If passed, plot to existing figure
    figsize: tuple (x, y)
        Figure size (if figure is None)

    Returns
    -------
    figure: PyVista figure
        Contains the plotted mesh

    """
    import pyvista as pv

    if figure is None:
        figure = pv.Plotter(window_size=figsize)
        figure.background_color = "white"
        pv.global_theme.font.color = "black"

    meshviz = pv.PolyData(
        mesh.vertices, np.hstack((np.repeat(3, len(mesh.faces))[:, None], mesh.faces))
    )
    figure.add_mesh(meshviz, **kwargs)
    figure.show(interactive_update=True)

    prop = pv.Property()
    if cull_back:
        prop.culling = "back"
    elif cull_front:
        prop.culling = "front"
    return figure


def polyline_from_points(points):
    import pyvista as pv

    poly = pv.PolyData()
    poly.points = points
    the_cell = np.arange(0, len(points), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(points))
    poly.lines = the_cell
    return poly


def plot_3d_current_loops(
    current_loops,
    colors="auto",
    figure=None,
    figsize=(800, 800),
    tube_radius=0.05,
    origin=np.array([0, 0, 0]),
):
    """
    Plot current loops (e.g. contour_polys given by scalar_contour()) in 3D using pyvista.

    Parameters
    ----------
    current_loops: list of (N_verts, 3) arrays with length N_loops
        List of segmented lines specifying the path of a current loop.
        Last vertex is assumed to be connected to the first.

    colors: 'auto' (default) or None or (3, ) tuple or  list of (3, ) tuples with same length as current_loops
        Optional specification of colors used for current direction. If auto (default),
        color blue-red based on current direction. If None, use middle-grey for all loops.
        If tuple, color all loops according to tuple. If list of tuples, color each loop accordingly.
    origin: array-like with length 3
        Shifts the origin used to determine the 'auto' coloring.
    tune_radius: float
        radius of loops plotted
    figure: existing pyvista figure
        Optional, if passed will plot to existing figure
    figsize: (x, y) tuple
        Optional, if plotting to new figure specifies the size (in pixels)

    Returns
    -------
    fig: pyvista figure

    """
    import pyvista as pv

    if figure is None:
        figure = pv.Plotter(window_size=figsize)
        figure.background_color = "white"
        pv.global_theme.font.color = "black"

    if colors is None:
        colors = "k" * len(current_loops)

    elif colors == "auto":
        colors = []

        palette = ["r", "b"]
        for loop_idx, loop in enumerate(current_loops):

            # Compute each loop segment
            segments = np.vstack(
                (loop[1:, :] - loop[0:-1, :], loop[0, :] - loop[-1, :])
            )

            # Find mean normal vector following right-hand rule, in loop centre
            centre_normal = np.mean(np.cross(segments, loop), axis=0)
            centre_normal /= np.linalg.norm(centre_normal, axis=-1)

            # Check if normal "points in" or "out" (towards or away from origin)
            origin_vector = np.mean(loop, axis=0) - origin

            colors.append(
                palette[int((np.sign(centre_normal @ origin_vector) + 1) / 2)]
            )

    elif isinstance(colors, str):
        colors = [colors] * len(current_loops)

    elif isinstance(colors, list):
        assert len(colors) == len(current_loops)

    else:
        raise ValueError("Invalid parameter for colors")

    for loop_idx, loop in enumerate(current_loops):

        polyline = polyline_from_points(loop[list(range(len(loop))) + [0]])

        tube = polyline.tube(radius=tube_radius)
        figure.add_mesh(tube, smooth_shading=True, color=colors[loop_idx])

        # Put two arrows on loop

        # First arrow on longest segment of loop
        longest_idx = np.argmax(
            np.linalg.norm(
                np.vstack((loop[1:, :] - loop[0:-1, :], loop[0, :] - loop[-1, :])),
                axis=-1,
            )
        )

        if longest_idx == len(loop) - 1:

            cone = pv.Cone(
                center=loop[-1, :].T,
                direction=(loop[0, :] - loop[-1, :]).T,
                height=0.5 * tube_radius / 0.05,
                radius=0.5 * tube_radius / 0.05 / 3,
                resolution=12,
            )

            figure.add_mesh(cone, color=colors[loop_idx])

        else:
            cone = pv.Cone(
                center=loop[longest_idx + 1, :].T,
                direction=(loop[longest_idx + 1, :] - loop[longest_idx, :]).T,
                height=0.5 * tube_radius / 0.05,
                radius=0.5 * tube_radius / 0.05 / 3,
                resolution=12,
            )

            figure.add_mesh(cone, color=colors[loop_idx])

    figure.view_isometric()
    figure.show(interactive_update=True)

    return figure


def plot_cross_section(
    X,
    Y,
    data,
    axes=None,
    cmap=None,
    colorbar=False,
    contours=10,
    log=False,
    vmin=None,
    vmax=None,
):
    """
    Plot scalar data on a plane

    Parameters
    ----------
    X: (N_x, N_y)
        meshgrid-generated 2D array with X coordinates for each data point
    Y: (N_x, N_y) array
        meshgrid-generated 2D array with Y coordinates for each data point
    data: (N_x, N_y) array
        2D array with data for each point in X, Y
    axes: matplotlib axes
        Optional, if passed plot to existing axes. If None (default), plot to new figure/axes
    cmap: str
        name of colormap to use for scalar data. If None (default), use viridis for
        all-positive/-negative data and RdBu otherwise
    colorbar: Boolean
        If True, plot colorbar for data
    contours: None, int or array-like
        If None, no contour lines are plotted. If int, plots specific number of
        equispaced contour lines. If array-like, plot contours at values in array.
    log: Boolean
        if True, colormap is log-scaled. Else, colormap is linear.
    vmin, vmax: float or None
        Explicit colormap minimum and maximum values. If not passed, range is according to data

    Returns
    -------
    axes: matplotlib axes with plot
    """

    if axes is None:
        fig, axes = plt.subplots(1, 1)

    # If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
    if cmap is None:
        if np.all(data > 0) or np.all(data < 0):
            cmap = "viridis"
        else:
            cmap = "RdBu"

    if vmin is None:
        vmin = np.min(data)

    if vmax is None:
        vmax = np.max(data)

    if log:
        norm = c.LogNorm(vmin, vmax)
    else:
        norm = c.Normalize(vmin, vmax)

    cont = axes.pcolormesh(X, Y, data, cmap=cmap, norm=norm, shading="gouraud")

    if colorbar:
        plt.colorbar(cont)

    if contours:
        clines = axes.contour(
            X,
            Y,
            data,
            levels=contours,
            norm=norm,
            antialiased=True,
            colors=("k",),
            linewidths=(1,),
        )

        axes.clabel(clines, fmt="%2.2f", colors="w", fontsize=10)

    axes.set_xlabel("X")
    axes.set_ylabel("Y")

    axes.figure.tight_layout()

    return axes


def plot_data_on_faces(
    mesh,
    data,
    figure=None,
    figsize=(800, 800),
    colorbar=False,
    vmin=None,
    vmax=None,
    **kwargs
):
    """Plot any data determined on the faces of a mesh

    Parameters
    ----------
    mesh: Trimesh mesh object
    data: (N_verts, ) array
        Scalar data to plot on mesh faces
    figure: existing pyvista figure
        Optional, if passed will plot to existing figure
    figsize: (x, y) tuple
        Optional, if plotting to new figure specifies the size (in pixels)
    colormap: str
        name of colormap to use for scalar data. If None (default), use viridis for
        all-positive/-negative data and RdBu otherwise
    colorbar: Boolean
        If True, plot colorbar for scalar data

    Returns
    -------
    fig: pyvista figure
    """
    import pyvista as pv

    if figure is None:
        figure = pv.Plotter(window_size=figsize)
        figure.background_color = "white"
        pv.global_theme.font.color = "black"

    # If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
    if "colormap" not in kwargs:
        if np.all(data > 0) or np.all(data < 0):
            kwargs["colormap"] = "viridis"
        else:
            kwargs["colormap"] = "RdBu"

    meshviz = pv.PolyData(
        mesh.vertices, np.hstack((np.repeat(3, len(mesh.faces))[:, None], mesh.faces))
    )

    if kwargs["colormap"] == "RdBu":
        rangemax = np.max(np.abs(data))
        if vmin is None:
            vmin = -rangemax * 1.01
        if vmax is None:
            vmax = rangemax * 1.01

    figure.add_mesh(
        meshviz, scalars=data, clim=[vmin, vmax], show_scalar_bar=colorbar, **kwargs
    )

    figure.show(interactive_update=True)

    return figure


def plot_data_on_vertices(
    mesh,
    data,
    figure=None,
    figsize=(800, 800),
    colorbar=False,
    ncolors=256,
    interpolate=True,
    cull_front=False,
    cull_back=False,
    autoscale=False,
    vmin=None,
    vmax=None,
    **kwargs
):
    """
    Plot scalar data defined on the vertices of a mesh.

    Parameters
    ----------
    mesh: Trimesh mesh object
    data: (N_verts, ) array
        Scalar data to plot on mesh vertices
    figure: existing pyvista figure
        Optional, if passed will plot to existing figure
    figsize: (x, y) tuple
        Optional, if plotting to new figure specifies the size (in pixels)
    cmap: str
        name of colormap to use for scalar data. If None (default), use viridis for
        all-positive/-negative data and RdBu otherwise
    colorbar: Boolean
        If True, plot colorbar for scalar data
    ncolors: int
        Number of colors to use
    interpolate: Boolean
        If True, interpolate scalar data for smoother look


    Returns
    -------
    fig: pyvista figure

    """

    import pyvista as pv

    if figure is None:
        figure = pv.Plotter(window_size=figsize)
        figure.background_color = "white"
        pv.global_theme.font.color = "black"

    # If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
    if "colormap" not in kwargs:
        if np.all(data > 0) or np.all(data < 0):
            kwargs["colormap"] = "viridis"
        else:
            kwargs["colormap"] = "RdBu"

    meshviz = pv.PolyData(
        mesh.vertices, np.hstack((np.repeat(3, len(mesh.faces))[:, None], mesh.faces))
    )

    if kwargs["colormap"] == "RdBu":
        rangemax = np.max(np.abs(data))
        if vmin is None:
            vmin = -rangemax * 1.01
        if vmax is None:
            vmax = rangemax * 1.01

    figure.add_mesh(
        meshviz,
        scalars=data,
        clim=[vmin, vmax],
        show_scalar_bar=colorbar,
        interpolate_before_map=interpolate,
        **kwargs
    )

    prop = pv.Property()
    if cull_back:
        prop.culling = "back"
    elif cull_front:
        prop.culling = "front"

    figure.show(interactive_update=True)

    return figure
