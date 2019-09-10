'''
Visualization functions tailored for bfieldtools.
Mainly wrappers and convenience helpers around mayavi and matplotlib functions
'''

from mayavi import mlab
import matplotlib.pyplot as plt
from matplotlib import colors

import numpy as np

def plot_3d_current_loops(current_loops, current_direction, colors=None, figure=None, figsize=(800, 800)):
    '''
    Plot current loops (e.g. contour_polys given by scalar_contour()) in 3D using mayavi.

    Parameters
    ----------
    current_loops: list of (N_verts, 3) arrays with length N_loops
        List of segmented lines specifying the path of a current loop.
        Last vertex is assumed to be connected to the first.
    current_direction: list of floats with length N_loops
        For each line/loop, specify the current direction by the sign of a scalar value.
    colors: list of (3, ) tuples with length 2
        Optional specification of colors used for current direction
    figure: existing mlab figure
        Optional, if passed will plot to existing figure
    figsize: (x, y) tuple
        Optional, if plotting to new figure specifies the size (in pixels)
    Returns
    -------
    fig: mlab figure

    '''

    if figure is None:
        fig = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=figsize)

    if colors is None:
        colors = [(1, 0, 0), (0, 0, 1)]

    for loop_idx, loop in enumerate(current_loops):
        mlab.plot3d(*loop[list(range(len(loop))) + [0]].T,
                    color=colors[int((np.sign(current_direction[loop_idx])+1)/2)],
                    tube_radius=0.05)

        mlab.quiver3d(*loop[1,:].T,
                  *(loop[1,:].T - loop[0,:].T),
                  mode='cone', scale_mode='none',
                  scale_factor=0.5,
                  color=colors[int((np.sign(current_direction[loop_idx])+1)/2)])

    fig.scene.isometric_view()

    return fig






def plot_cross_section(X, Y, data, axes=None, cmap=None, colorbar=False, contours=10, log=False, vmin=None, vmax=None):
    '''
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
    '''


    if axes is None:
        fig, axes = plt.subplots(1,1)

    #If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
    if cmap is None:
        if np.all(data > 0) or np.all(data < 0):
            cmap='viridis'
        else:
            cmap='RdBu'

    if log:
        norm = colors.LogNorm()
    else:
        norm = colors.Normalize()

    if vmin is None:
        vmin = np.min(data)

    if vmax is None:
        vmax = np.max(data)

    cont = axes.pcolormesh(X, Y,
                         data,
                         cmap=cmap,
                         vmin=vmin,
                         vmax=vmax,
                         norm=norm,
                         shading='gouraud')

    if contours:
        clines = axes.contour(X, Y,
                            data,
                            levels=contours,
                            norm=norm,
                            antialiased=True,
                            colors=('k',),
                            linewidths=(1,))

        axes.clabel(clines, fmt='%2.2f', colors='w', fontsize=10)

    axes.set_xlabel('X')
    axes.set_ylabel('Y')

    axes.figure.tight_layout()

    return axes


#def plot_field_falloff(axis, points, mesh, current_density, figsize):
#    '''
#
#
#    '''
#
#    if axes is None:
#        fig, axes = plt.subplots(1,1)
#
#    #If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
#    if cmap is None:
#        if np.all(data > 0) or np.all(data < 0):
#            cmap='viridis'
#        else:
#            cmap='RdBu'
#
#    r = np.zeros(())
#
#    compute_C
#
#    fig = plt.figure(figsize=figsize)
#
#    plt.semilogy(label=)
#
#
#    return fig

def plot_data_on_faces(mesh, data, figure=None, figsize=(800, 800), cmap=None, colorbar=False, ncolors=32):
    ''' Plot any data determined on the faces of a mesh

        Parameters
        ----------
        mesh: Trimesh mesh object
        data: (N_verts, ) array
            Scalar data to plot on mesh faces
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

        Returns
        -------
        fig: mlab figure
    '''


    if figure is None:
        fig = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=figsize)

    #If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
    if cmap is None:
        if np.all(data > 0) or np.all(data < 0):
            cmap='viridis'
        else:
            cmap='RdBu'

    v = mesh.vertices
    f = mesh.faces

    s = mlab.pipeline.triangular_mesh_source(*v.T, f)
    s.mlab_source.dataset.cell_data.scalars = data

    s.mlab_source.dataset.cell_data.scalars.name = 'Cell data'

    s.mlab_source.update()
    s2 = mlab.pipeline.set_active_attribute(s,cell_scalars='Cell data')
    surf = mlab.pipeline.surface(s2)

    lutmanager = surf.parent.scalar_lut_manager
    lutmanager.lut_mode = cmap

    if cmap == 'RdBu':
        rangemax = np.max(np.abs(data))
        lutmanager.data_range = np.array([-rangemax*1.01,rangemax*1.01])
        lutmanager.use_default_range = False

    lutmanager.number_of_colors = ncolors

    return fig


def plot_data_on_vertices(mesh, data, figure=None, figsize=(800, 800), cmap=None, colorbar=False, ncolors=32, interpolate=True):
    '''
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
    Returns
    -------
    fig: mlab figure

    '''

    if figure is None:
        fig = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=figsize)

    #If data is all-positive or all-negative, use viridis. Otherwise, use Red-Blue colormap
    if cmap is None:
        if np.all(data > 0) or np.all(data < 0):
            cmap='viridis'
        else:
            cmap='RdBu'


    surf = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=data, colormap=cmap)
    mlab.colorbar(surf)
    surf.actor.mapper.interpolate_scalars_before_mapping = interpolate

    lutmanager = surf.parent.scalar_lut_manager
    lutmanager.number_of_colors = ncolors

    if cmap == 'RdBu':
        rangemax = np.max(np.abs(data))
        lutmanager.data_range = np.array([-rangemax*1.01,rangemax*1.01])
        lutmanager.use_default_range = False

    return fig


