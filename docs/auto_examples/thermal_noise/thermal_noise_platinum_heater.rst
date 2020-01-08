.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_thermal_noise_thermal_noise_platinum_heater.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_thermal_noise_thermal_noise_platinum_heater.py:


Thin-film heater noise
=========================

This example computes the thermal magnetic noise produced by a platinum
thin-film heater geometry.


.. code-block:: default


    import numpy as np
    import matplotlib.pyplot as plt
    import trimesh
    from mayavi import mlab

    from bfieldtools.thermal_noise import compute_current_modes, compute_dc_Bnoise

    import pkg_resources








Fix the simulation parameters and load the heater geometry


.. code-block:: default


    d = 500e-9 #Film thickness in meters
    sigma = 1/(16.592 * 1e-6 * 1e-2) #Platinum @ 450 K
    T = 450 #Temperature in K
    kB = 1.38064852e-23 #Boltzman constant
    mu0 = 4*np.pi*1e-7


    scale_factor = 1e-3



    mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/thin_film_heater.stl'))

    #Subdivide mesh for higher accuracy if needed
    #mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    # Center the mesh at the origin, apply scaling
    mesh.apply_scale(scale_factor)

    mesh.vertices[:, 2] = 0
    mesh.vertices[:, 1] -= np.mean(mesh.vertices[:, 1])
    mesh.vertices[:, 0] -= np.mean(mesh.vertices[:, 0])







Visualize the geometry.


.. code-block:: default


    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(800, 800))

    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)




.. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_platinum_heater_001.png
    :class: sphx-glr-single-img




Compute the normalized thermal current modes, and thereafter compute the
 magnetic field noise caused by the currents. Finally, visualize the result.


.. code-block:: default


    vl = compute_current_modes(mesh)

    Np = 30

    zl = np.linspace(0.1, 5, Np) * scale_factor
    fp = np.array((np.zeros(zl.shape), np.zeros(zl.shape)-0.001, zl)).T

    B = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)

    fig = plt.figure(figsize=(6, 4))

    plt.semilogy(zl*1e3, np.linalg.norm(B, axis=1)*1e15, 'k')
    plt.xlabel('Distance (mm)')
    plt.ylabel('DC noise amplitude (fT/rHz)')

    plt.grid()
    plt.title('Thermal noise falloff')
    fig.tight_layout()





.. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_platinum_heater_002.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 1011 vertices by 30 target points... took 0.03 seconds.



Compute the field on a 3D grid and visualize isosurfaces.


.. code-block:: default


    plane_extent = 3.5
    Ngrid = 40

    xx = np.linspace(-plane_extent, plane_extent, Ngrid) * scale_factor
    yy = np.linspace(-plane_extent, plane_extent, Ngrid) * scale_factor
    zz = np.array([0.1, 0.25, 0.5, 1, 1.5]) * scale_factor
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    grid_points = np.vstack((x, y, z)).T




    B_grid = compute_dc_Bnoise(mesh, vl, grid_points, sigma, d, T)

    B_grid_matrix = B_grid.reshape((Ngrid, Ngrid, len(zz), 3))

    B_grid_matrix_norm = np.linalg.norm(B_grid_matrix, axis=-1)


    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(800, 800))

    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)

    field = mlab.pipeline.vector_field(X, Y, Z, B_grid_matrix[:,:,:,0], B_grid_matrix[:,:,:,1], B_grid_matrix[:,:,:,2],
                                  scalars=B_grid_matrix_norm, name='B-field')


    iso = mlab.pipeline.iso_surface(field,
                                    opacity=0.3,
                                    colormap='viridis',
                                    contours=[20e-15, 5e-15, 1e-15, 1e-16],
                                    vmax=20e-15,
                                    vmin=1e-16)

    # A trick to make transparency look better: cull the front face
    iso.actor.property.frontface_culling = False




.. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_platinum_heater_003.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 1011 vertices by 8000 target points... took 1.53 seconds.



Plot the noise level at horizontal planes at different distance.


.. code-block:: default




    from matplotlib import colors

    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(20, 4))
    axes = axes.flatten()
    B_scale = 1e15

    for ax_idx, ax in enumerate(axes):
        cont = ax.pcolormesh(X[:, :, ax_idx]*scale_factor, Y[:, :, ax_idx]*scale_factor,
                             B_scale * B_grid_matrix_norm[:, :, ax_idx],
                             cmap='viridis',
                             vmin=B_scale * 1e-17, vmax=B_scale * 5e-14,
                             norm=colors.LogNorm(),
                             shading='gouraud')

        clines = ax.contour(X[:, :, ax_idx]*scale_factor, Y[:, :, ax_idx]*scale_factor, B_scale * B_grid_matrix_norm[:, :, ax_idx],
                            levels=B_scale * np.array([1e-17, 5e-17, 1e-16, 5e-16, 1e-15, 2.5e-15, 5e-15, 1e-14, 2.5e-14, 5e-14]),
                            norm=colors.LogNorm(),
                            antialiased=True,
                            colors=('k',),
                            linewidths=(3,))
        ax.clabel(clines, fmt='%2.2f', colors='w', fontsize=10)

        ax.set_title('Distance %.2f mm'%(Z[0, 0, ax_idx]*1e3))
        ax.set_xlabel('(mm)')
        ax.set_ylabel('(mm)')

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    fig.tight_layout()

    fig.subplots_adjust(right=0.925)
    cbar_ax = fig.add_axes([0.95, 0.15, 0.01, 0.7])
    cbar = fig.colorbar(cont, cax=cbar_ax)
    cbar.set_label('DC magnetic field noise amplitude (fT/rHz)')


.. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_platinum_heater_004.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  19.595 seconds)

**Estimated memory usage:**  962 MB


.. _sphx_glr_download_auto_examples_thermal_noise_thermal_noise_platinum_heater.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: thermal_noise_platinum_heater.py <thermal_noise_platinum_heater.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: thermal_noise_platinum_heater.ipynb <thermal_noise_platinum_heater.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
