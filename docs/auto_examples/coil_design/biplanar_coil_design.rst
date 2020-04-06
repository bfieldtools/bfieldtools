.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_coil_design_biplanar_coil_design.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_coil_design_biplanar_coil_design.py:


Biplanar coil design
====================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes.


.. code-block:: default


    import numpy as np
    import matplotlib.pyplot as plt
    from mayavi import mlab
    import trimesh


    from bfieldtools.conductor import Conductor
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

    planemesh.apply_scale(scaling_factor*1.6)

    #Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 5, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                             planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff, 
                         planemesh.faces, process=False)

    joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = Conductor(mesh_obj=joined_planes, fix_normals=True, basis_name='suh', N_suh=100)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Calculating surface harmonics expansion...
    Computing the laplacian matrix...
    Computing the mass matrix...




Set up target and stray field points


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 0]) * scaling_factor

    sidelength = 2 * scaling_factor
    n = 8
    xx = np.linspace(-sidelength/2, sidelength/2, n)
    yy = np.linspace(-sidelength/2, sidelength/2, n)
    zz = np.linspace(-sidelength/2, sidelength/2, n)
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    target_points = np.array([x, y, z]).T

    #Turn cube into sphere by rejecting points "in the corners"
    target_points = target_points[np.linalg.norm(target_points, axis=1) < sidelength/2]  + center



    #    #Here, the stray field points are on a spherical surface
    stray_radius = 20 * scaling_factor
    #    stray_length = 20 * scaling_factor
    #
    #    stray_points = cylinder_points(radius=stray_radius,
    #                                   length = stray_length,
    #                                   nlength = 5,
    #                                   nalpha = 30,
    #                                   orientation=np.array([1, 0, 0]))
    #
    stray_points_mesh = trimesh.creation.icosphere(subdivisions=3, radius=stray_radius)
    stray_points = stray_points_mesh.vertices + center

    n_stray_points = len(stray_points)









Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function

    target_field = np.zeros(target_points.shape)
    target_field[:, 0] = target_field[:, 0] + 1

    target_spec = {'coupling':coil.B_coupling(target_points), 'abs_error':0.01, 'target':target_field}
    stray_spec = {'coupling':coil.B_coupling(stray_points), 'abs_error':0.001, 'target':np.zeros((n_stray_points, 3))}

    bfield_specification = [target_spec, stray_spec]





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.32 seconds.
    Computing magnetic field coupling matrix, 3184 vertices by 642 target points... took 0.92 seconds.




# Compute the optimal stream function, either using a numerical solver or regularized least squares


.. code-block:: default

    import mosek

    coil.s, prob = optimize_streamfunctions(coil,
                                       [target_spec, stray_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 120 chunks (5898 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 61.14 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 1062            
      Cones                  : 1               
      Scalar variables       : 203             
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 1062            
      Cones                  : 1               
      Scalar variables       : 203             
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 101
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 1062              conic                  : 102             
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.02              dense det. time        : 0.00            
    Factor     - ML order time          : 0.00              GP order time          : 0.00            
    Factor     - nonzeros before factor : 5151              after factor           : 5151            
    Factor     - dense dim.             : 0                 flops                  : 6.22e+06        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.6e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  0.05  
    1   9.4e+00  5.8e-01  1.0e+00  -3.70e-01  5.663010048e+00   5.041518979e+00   5.8e-01  0.05  
    2   4.1e+00  2.5e-01  9.3e-02  -1.70e-01  3.145612954e+01   3.097516480e+01   2.5e-01  0.05  
    3   3.3e+00  2.0e-01  6.4e-02  7.92e-01   4.958926055e+01   4.927719670e+01   2.0e-01  0.05  
    4   1.3e+00  8.3e-02  2.6e-02  1.97e+00   6.930700588e+01   6.922086898e+01   8.3e-02  0.05  
    5   9.3e-01  5.7e-02  1.6e-02  1.06e+00   7.579364951e+01   7.573414764e+01   5.7e-02  0.05  
    6   4.8e-02  3.0e-03  1.7e-04  1.16e+00   8.585002834e+01   8.584707817e+01   3.0e-03  0.06  
    7   4.2e-03  2.6e-04  5.3e-06  1.02e+00   8.636391171e+01   8.636369220e+01   2.6e-04  0.06  
    8   2.3e-07  1.4e-08  1.9e-12  1.00e+00   8.641835039e+01   8.641835038e+01   1.4e-08  0.06  
    Optimizer terminated. Time: 0.06    


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 8.6418350393e+01    nrm: 2e+02    Viol.  con: 9e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 8.6418350377e+01    nrm: 2e+02    Viol.  con: 3e-08    var: 1e-10    cones: 0e+00  




Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.s, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.s

    mlab.quiver3d(*target_points.T, *B_target.T)








.. image:: /auto_examples/coil_design/images/sphx_glr_biplanar_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.vectors.Vectors object at 0x00000274B81181A8>



Plot cross-section of magnetic field and magnetic potential of the discretized loops


.. code-block:: default


    from bfieldtools.line_magnetics import magnetic_field, scalar_potential

    x = y = np.linspace(-12, 12, 250)
    X,Y = np.meshgrid(x, y, indexing='ij')
    points = np.zeros((X.flatten().shape[0], 3))
    points[:, 0] = X.flatten()
    points[:, 1] = Y.flatten()

    B = np.zeros_like(points)
    U = np.zeros((points.shape[0],))
    for loop_idx in range(len(loops)):
        B += magnetic_field(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)
        U += scalar_potential(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)

    B = B.T[:2].reshape(2, x.shape[0], y.shape[0])
    lw = np.sqrt(B[0]**2 + B[1]**2)
    lw = 2*lw/np.max(lw)

    U = U.reshape(x.shape[0], y.shape[0])

    plt.figure()

    plt.pcolormesh(X, Y, U.T, cmap='seismic', shading='gouraud')
    #plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
    #           extent=(x.min(), x.max(), y.min(), y.max()))

    seed_points=points[:,:2]*0.3

    plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k', integration_direction='both',
                   start_points=seed_points)
    plt.axis('equal')
    plt.axis('off')
    for loop in loops:
        plt.plot(loop[:,1], loop[:,0], 'k', linewidth=4, alpha=0.1)

    plt.tight_layout()





.. image:: /auto_examples/coil_design/images/sphx_glr_biplanar_coil_design_002.png
    :class: sphx-glr-single-img





Lets also do the same coil optimization using regularized least-squares


.. code-block:: default



    from bfieldtools.coil_optimize import optimize_lsq

    coil.s2 = optimize_lsq(coil,
                          [target_spec, stray_spec],
                          objective='minimum_resistive_energy',
                          reg=1e6)








.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing the resistance matrix...
    Error tolerances in specification will be ignored when use lsq




Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.s2, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.s2

    mlab.quiver3d(*target_points.T, *B_target.T)








.. image:: /auto_examples/coil_design/images/sphx_glr_biplanar_coil_design_003.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.vectors.Vectors object at 0x00000274B83E36D0>



Plot cross-section of magnetic field and magnetic potential of the discretized loops


.. code-block:: default


    from bfieldtools.line_magnetics import magnetic_field, scalar_potential

    x = y = np.linspace(-12, 12, 250)
    X,Y = np.meshgrid(x, y, indexing='ij')
    points = np.zeros((X.flatten().shape[0], 3))
    points[:, 0] = X.flatten()
    points[:, 1] = Y.flatten()

    B = np.zeros_like(points)
    U = np.zeros((points.shape[0],))
    for loop_idx in range(len(loops)):
        B += magnetic_field(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)
        U += scalar_potential(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)

    B = B.T[:2].reshape(2, x.shape[0], y.shape[0])
    lw = np.sqrt(B[0]**2 + B[1]**2)
    lw = 2*lw/np.max(lw)

    U = U.reshape(x.shape[0], y.shape[0])

    plt.figure()

    plt.pcolormesh(X, Y, U.T, cmap='seismic', shading='gouraud')
    #plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
    #           extent=(x.min(), x.max(), y.min(), y.max()))

    seed_points=points[:,:2]*0.3

    plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k', integration_direction='both',
                   start_points=seed_points)
    plt.axis('equal')
    plt.axis('off')
    for loop in loops:
        plt.plot(loop[:,1], loop[:,0], 'k', linewidth=4, alpha=0.1)

    plt.tight_layout()




.. image:: /auto_examples/coil_design/images/sphx_glr_biplanar_coil_design_004.png
    :class: sphx-glr-single-img






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 4 minutes  14.097 seconds)


.. _sphx_glr_download_auto_examples_coil_design_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: biplanar_coil_design.py <biplanar_coil_design.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: biplanar_coil_design.ipynb <biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
