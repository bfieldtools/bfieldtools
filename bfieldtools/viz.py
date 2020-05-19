"""
Visualization functions tailored for bfieldtools.
Mainly wrappers and convenience helpers around mayavi and matplotlib functions
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
    Plot the mesh surface in mayavi.

    Parameters
    ----------
    mesh: Trimesh mesh object
        mesh to be plotted
    cull_front: Boolean (False)
        If True, cull front of mesh
    cull_back: Boolean (False)
        If True, cull back of mesh
    figure: mayavi figure or None
        If passed, plot to existing figure
    figsize: tuple (x, y)
        Figure size (if figure is None)

    Returns
    -------
    figure: mayavi figure
        Contains the plotted mesh

    """
    from mayavi import mlab

    if figure is None:
        figure = mlab.figure(
            None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=figsize
        )

    meshviz = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, **kwargs)

    meshviz.actor.mapper.scalar_visibility = False

    meshviz.actor.property.frontface_culling = cull_front
    meshviz.actor.property.backface_culling = cull_back
    return figure


def plot_3d_current_loops(
    current_loops,
    colors="auto",
    figure=None,
    figsize=(800, 800),
    tube_radius=0.05,
    origin=np.array([0, 0, 0]),
):
    """
    Plot current loops (e.g. contour_polys given by scalar_contour()) in 3D using mayavi.

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
    figure: existing mlab figure
        Optional, if passed will plot to existing figure
    figsize: (x, y) tuple
        Optional, if plotting to new figure specifies the size (in pixels)

    Returns
    -------
    fig: mlab figure

    """
    from mayavi import mlab

    if figure is None:
        figure = mlab.figure(
            None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=figsize
        )

    if colors is None:
        colors = [(0.5, 0.5, 0.5)] * len(current_loops)

    elif isinstance(colors, tuple):
        colors = [colors] * len(current_loops)

    elif isinstance(colors, list):
        assert len(colors) == len(current_loops)

    elif colors == "auto":
        colors = []

        palette = [(1, 0, 0), (0, 0, 1)]
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
    else:
        raise ValueError("Invalid parameter for colors")

    for loop_idx, loop in enumerate(current_loops):
        mlab.plot3d(
            *loop[list(range(len(loop))) + [0]].T,
            color=colors[loop_idx],
            tube_radius=tube_radius
        )

        # Put two arrows on loop

        # First arrow on longest segment of loop
        longest_idx = np.argmax(
            np.linalg.norm(
                np.vstack((loop[1:, :] - loop[0:-1, :], loop[0, :] - loop[-1, :])),
                axis=-1,
            )
        )

        if longest_idx == len(loop) - 1:
            arrow1 = mlab.quiver3d(
                *loop[-1, :].T,
                *(loop[0, :] - loop[-1, :]).T,
                mode="cone",
                scale_mode="none",
                scale_factor=0.5 * tube_radius / 0.05,
                color=colors[loop_idx]
            )
        else:
            arrow1 = mlab.quiver3d(
                *loop[longest_idx + 1, :].T,
                *(loop[longest_idx + 1, :] - loop[longest_idx, :]).T,
                mode="cone",
                scale_mode="none",
                scale_factor=0.5 * tube_radius / 0.05,
                color=colors[loop_idx]
            )

        arrow1.glyph.glyph_source.glyph_position = "center"
        arrow1.glyph.glyph_source.glyph_source.radius = 0.3
        arrow1.glyph.glyph_source.glyph_source.height = 0.5

    #        #Second arrow on the element "half away"
    #
    #        opposite_idx = int((longest_idx + len(loop)/2) % len(loop))
    #
    #        if opposite_idx == len(loop)-1:
    #            arrow1 = mlab.quiver3d(*loop[-1,:].T,
    #                      *(loop[0,:] - loop[-1,:]).T,
    #                      mode='cone', scale_mode='none',
    #                      scale_factor=0.5 * tube_radius/0.05,
    #                      color=colors[loop_idx])
    #        else:
    #            arrow2 = mlab.quiver3d(*loop[opposite_idx+1,:].T,
    #                      *(loop[opposite_idx+1,:] - loop[opposite_idx,:]).T,
    #                      mode='cone', scale_mode='none',
    #                      scale_factor=0.5 * tube_radius/0.05,
    #                      color=colors[loop_idx])
    #        arrow2.glyph.glyph_source.glyph_position = 'center'
    #        arrow2.glyph.glyph_source.glyph_source.radius = 0.3
    #        arrow2.glyph.glyph_source.glyph_source.height = 0.5

    figure.scene.isometric_view()

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

    if log:
        norm = c.LogNorm()
    else:
        norm = c.Normalize()

    if vmin is None:
        vmin = np.min(data)

    if vmax is None:
        vmax = np.max(data)

    cont = axes.pcolormesh(
        X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm, shading="gouraud"
    )

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
    ncolors=256,
    vmin=None,
    vmax=None,
    **kwargs
):
    """ Plot any data determined on the faces of a mesh

        Parameters
        ----------
        mesh: Trimesh mesh object
        data: (N_verts, ) array
            Scalar data to plot on mesh faces
        figure: existing mlab figure
            Optional, if passed will plot to existing figure
        figsize: (x, y) tuple
            Optional, if plotting to new figure specifies the size (in pixels)
        colormap: str
            name of colormap to use for scalar data. If None (default), use viridis for
            all-positive/-negative data and RdBu otherwise
        colorbar: Boolean
            If True, plot colorbar for scalar data
        ncolors: int
            Number of colors to use

        Returns
        -------
        fig: mlab figure
    """
    from mayavi import mlab

    if figure is None:
        figure = mlab.figure(
            None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=figsize
        )

    # If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
    if "colormap" not in kwargs:
        if np.all(data > 0) or np.all(data < 0):
            kwargs["colormap"] = "viridis"
        else:
            kwargs["colormap"] = "RdBu"

    v = mesh.vertices
    f = mesh.faces

    s = mlab.pipeline.triangular_mesh_source(*v.T, f)
    s.mlab_source.dataset.cell_data.scalars = data

    s.mlab_source.dataset.cell_data.scalars.name = "Cell data"

    s.mlab_source.update()
    s2 = mlab.pipeline.set_active_attribute(s, cell_scalars="Cell data")
    surf = mlab.pipeline.surface(s2, **kwargs)

    if colorbar:
        mlab.colorbar(surf)

    lutmanager = surf.parent.scalar_lut_manager
    lutmanager.lut_mode = kwargs["colormap"]

    if kwargs["colormap"] == "RdBu":
        rangemax = np.max(np.abs(data))
        if vmin is None:
            vmin = -rangemax * 1.01
        if vmax is None:
            vmax = rangemax * 1.01
        lutmanager.data_range = np.array([vmin, vmax])
        lutmanager.use_default_range = False

    lutmanager.number_of_colors = ncolors

    return surf


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
    figure: existing mlab figure
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
    opacity: float
        Opacity of rendered mesh

    Returns
    -------
    fig: mlab figure

    """
    from mayavi import mlab

    if figure is None:
        figure = mlab.figure(
            None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=figsize
        )

    # If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
    if "colormap" not in kwargs:
        if np.all(data > 0) or np.all(data < 0):
            kwargs["colormap"] = "viridis"
        else:
            kwargs["colormap"] = "RdBu"

    surf = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=data, **kwargs)
    surf.actor.property.frontface_culling = cull_front
    surf.actor.property.backface_culling = cull_back
    if colorbar:
        mlab.colorbar(surf)
    surf.actor.mapper.interpolate_scalars_before_mapping = interpolate

    lutmanager = surf.module_manager.scalar_lut_manager
    lutmanager.number_of_colors = ncolors

    if (np.all(data > 0) or np.all(data < 0)) and autoscale:
        rangemax = np.max(np.abs(data))
        lutmanager.data_range = np.array([-rangemax * 1.01, rangemax * 1.01])
        lutmanager.use_default_range = False
    elif vmax is not None:
        if vmin is None:
            lutmanager.data_range = np.array([-vmax, vmax])
        else:
            lutmanager.data_range = np.array([vmin, vmax])
        lutmanager.use_default_range = False

    return surf
