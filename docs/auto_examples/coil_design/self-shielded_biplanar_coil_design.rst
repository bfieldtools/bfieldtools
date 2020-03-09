.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_self-shielded_biplanar_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_self-shielded_biplanar_coil_design.py:


Analytical self-shielded biplanar coil design
==============================================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes. In addition, the coils have an outer surface
for which (in a linear fashion) a secondary current is created, which zeroes the
normal component of the field produced by the primary coil at the secondary coil
surface. The combination of the primary and secondary coil currents are specified to create
the target field, and their combined inductive energy is minimized.

NB. The secondary coil current is entirely a function of the primary coil current
and the geometry.


.. code-block:: default


    import numpy as np
    import matplotlib.pyplot as plt
    from mayavi import mlab
    import trimesh

    from bfieldtools.mesh_class import Conductor, StreamFunction
    from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic, scalar_potential_coupling
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops, plot_data_on_vertices

    import pkg_resources
    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 0.1

    #Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

    planemesh.apply_scale(scaling_factor)

    #Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 4, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                             planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                         planemesh.faces, process=False)

    mesh1 = coil_plus.union(coil_minus)
    mesh2 = mesh1.copy()
    mesh2.apply_scale(1.4)

    coil = Conductor(mesh_obj=mesh1, basis_name = 'inner', N_sph = 4)
    shieldcoil = Conductor(mesh_obj=mesh2, basis_name = 'inner', N_sph = 4)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    face_normals didn't match triangles, ignoring!



Plot geometry


.. code-block:: default

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    coil.plot_mesh(opacity=0.2, figure=f)
    shieldcoil.plot_mesh(opacity=0.2, figure=f)




.. image:: /auto_examples/coil_design/images/sphx_glr_self-shielded_biplanar_coil_design_001.png
    :class: sphx-glr-single-img




Compute inductances and coupling


.. code-block:: default



    M11 = coil.inductance
    M22 = shieldcoil.inductance
    M21 = shieldcoil.mutual_inductance(coil)


    # Mapping from I1 to I2, constraining flux through shieldcoil to zero
    P = -np.linalg.solve(M22, M21)

    A1, Beta1 = coil.sph_couplings
    A2, Beta2 = shieldcoil.sph_couplings





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 80 chunks (11499 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 38.97 seconds.
    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 80 chunks (11289 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 39.28 seconds.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 80 chunks (11177 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Computing coupling matrices
    l = 1 computed
    l = 2 computed
    l = 3 computed
    l = 4 computed
    Computing coupling matrices
    l = 1 computed
    l = 2 computed
    l = 3 computed
    l = 4 computed



Precalculations for the solution


.. code-block:: default


    # Minimization of magnetic energy with spherical harmonic constraint
    C = Beta1 + Beta2 @ P
    M = M11 + M21.T @ P

    #Regularization
    from scipy.linalg import eigvalsh
    ssmax = eigvalsh(C.T @ C, M, eigvals=[M.shape[1]-1, M.shape[1]-1])







Specify spherical harmonic and calculate corresponding shielded field


.. code-block:: default

    beta = np.zeros(Beta1.shape[0])
    #beta[7] = 1 # Gradient
    beta[2] = 1 # Homogeneous

    # Minimum residual
    _lambda=1e3
    # Minimum energy
    #_lambda=1e-3
    I1inner = np.linalg.solve(C.T @ C + M*ssmax/_lambda, C.T @ beta)

    I2inner = P @ I1inner

    coil.s = StreamFunction(I1inner, coil)
    shieldcoil.s = StreamFunction(I2inner, shieldcoil)

    #s = mlab.triangular_mesh(*mesh1.vertices.T, mesh1.faces, scalars=I1)
    #s.enable_contours=True
    #s = mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=I2)
    #s.enable_contours=True









.. code-block:: default


    x = y = np.linspace(-0.8, 0.8, 150)
    X,Y = np.meshgrid(x, y, indexing='ij')
    points = np.zeros((X.flatten().shape[0], 3))
    points[:, 0] = X.flatten()
    points[:, 1] = Y.flatten()


    CB1 = coil.B_coupling(points)
    CB2 = shieldcoil.B_coupling(points)

    CU1 = coil.U_coupling(points)
    CU2 = shieldcoil.U_coupling(points)

    B1 = CB1 @ coil.s
    B2 = CB2 @ shieldcoil.s

    U1 = CU1 @ coil.s
    U2 = CU2 @ shieldcoil.s







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 22500 target points... took 24.04 seconds.
    Computing magnetic field coupling matrix, 3184 vertices by 22500 target points... took 23.25 seconds.
    Computing scalar potential coupling matrix, 3184 vertices by 22500 target points... took 106.70 seconds.
    Computing scalar potential coupling matrix, 3184 vertices by 22500 target points... took 106.34 seconds.



Now, plot the field streamlines and scalar potential


.. code-block:: default

    cc1 = scalar_contour(mesh1, mesh1.vertices[:,2], contours= [-0.001])[0]
    cc2 = scalar_contour(mesh2, mesh2.vertices[:,2], contours= [-0.001])[0]
    cx10 = cc1[0][:,1]
    cy10 = cc1[0][:,0]
    cx20 = cc2[0][:,1]
    cy20 = cc2[0][:,0]

    cx11 = np.vstack(cc1[1:])[:,1]
    cy11 = np.vstack(cc1[1:])[:,0]
    cx21 = np.vstack(cc2[1:])[:,1]
    cy21 = np.vstack(cc2[1:])[:,0]

    B = (B1.T + B2.T)[:2].reshape(2, x.shape[0], y.shape[0])
    lw = np.sqrt(B[0]**2 + B[1]**2)
    lw = 2*np.log(lw/np.max(lw)*np.e+1.1)

    xx = np.linspace(-1,1, 16)
    #seed_points = 0.56*np.array([xx, -np.sqrt(1-xx**2)])
    #seed_points = np.hstack([seed_points, (0.56*np.array([xx, np.sqrt(1-xx**2)]))])
    #seed_points = np.hstack([seed_points, (0.56*np.array([np.zeros_like(xx), xx]))])
    seed_points = np.array([cx10+0.001, cy10])
    seed_points = np.hstack([seed_points, np.array([cx11-0.001, cy11])])
    seed_points = np.hstack([seed_points, (0.56*np.array([np.zeros_like(xx), xx]))])

    #plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
    #               start_points=seed_points.T, integration_direction='both')
    U = (U1 + U2).reshape(x.shape[0], y.shape[0])
    U /= np.max(U)
    plt.figure()
    plt.contourf(X,Y, U.T, cmap='seismic', levels=40)
    #plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
    #           extent=(x.min(), x.max(), y.min(), y.max()))
    plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
                   start_points=seed_points.T, integration_direction='both',
                   arrowsize=0.1)

    #plt.plot(seed_points[0], seed_points[1], '*')

    plt.plot(cx10, cy10, linewidth=3.0, color='gray')
    plt.plot(cx20, cy20, linewidth=3.0, color='gray')
    plt.plot(cx11, cy11, linewidth=3.0, color='gray')
    plt.plot(cx21, cy21, linewidth=3.0, color='gray')
    plt.axis('image')

    plt.xticks([])
    plt.yticks([])





.. image:: /auto_examples/coil_design/images/sphx_glr_self-shielded_biplanar_coil_design_002.png
    :class: sphx-glr-single-img




Do a quick 3D plot


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))

    coil.s.plot(figure=f, contours=20)
    shieldcoil.s.plot(figure=f, contours=20)


.. image:: /auto_examples/coil_design/images/sphx_glr_self-shielded_biplanar_coil_design_003.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 11 minutes  54.460 seconds)

**Estimated memory usage:**  11658 MB


.. _sphx_glr_download_auto_examples_coil_design_self-shielded_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: self-shielded_biplanar_coil_design.py <self-shielded_biplanar_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: self-shielded_biplanar_coil_design.ipynb <self-shielded_biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
