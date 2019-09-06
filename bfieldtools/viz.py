'''
Visualization functions tailored for bfieldtools.
Mainly wrappers and convenience helpers around mayavi and matplotlib functions
'''

from mayavi import mlab
import matplotlib.pyplot as plt

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


def plot_scalar_on_mesh(mesh, scalar):
    '''

    '''
    return


def plot_cross_section(figsize):

    fig = plt.figure(figsize=figsize)
    return fig


def plot_field_falloff(axis, points, mesh, current_density, figsize):
    '''


    '''

    for

    r = np.zeros(())

    compute_C

    fig = plt.figure(figsize=figsize)

    plt.semilogy(label=)


    return fig


