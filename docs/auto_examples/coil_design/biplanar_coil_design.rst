.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_biplanar_coil_design.py>` to download the full example code
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

    from bfieldtools.mesh_class import MeshWrapper
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
    coil = MeshWrapper(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)







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

    target_rel_error = np.zeros_like(target_field)
    target_rel_error[:, 0] += 0.01

    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 0] += 0.001
    target_abs_error[:, 1:3] += 0.005

    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}
    stray_spec = {'coupling':coil.B_coupling(stray_points), 'abs_error':0.01, 'rel_error':0, 'target':np.zeros((n_stray_points, 3))}

    bfield_specification = [target_spec, stray_spec]





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.30 seconds.
    Computing magnetic field coupling matrix, 3184 vertices by 642 target points... took 0.72 seconds.



Run QP solver


.. code-block:: default

    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       [target_spec, stray_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 405514 MiB required for 3184 times 3184 vertices...
    Computing inductance matrix in 42 chunks since 9694 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 83.09 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 7710            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 7710            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 7710              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 1.84              dense det. time        : 0.00            
    Factor     - ML order time          : 0.15              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 5.26e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   4.1e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  93.50 
    1   2.5e+01  6.2e-01  1.4e+00  -7.89e-01  2.213390980e+00   1.706833683e+00   6.2e-01  94.09 
    2   1.4e+01  3.5e-01  9.0e-01  -6.22e-01  1.906120449e+01   1.915662410e+01   3.5e-01  94.66 
    3   9.3e+00  2.3e-01  5.9e-01  -3.63e-01  9.015555708e+01   9.053356335e+01   2.3e-01  95.21 
    4   5.2e+00  1.3e-01  3.2e-01  -1.44e-01  1.375282726e+02   1.381228711e+02   1.3e-01  95.75 
    5   2.4e+00  5.8e-02  1.2e-01  1.85e-01   3.300859071e+02   3.305769736e+02   5.8e-02  96.30 
    6   1.7e+00  4.2e-02  7.7e-02  5.49e-01   3.938765068e+02   3.942839359e+02   4.2e-02  96.85 
    7   5.0e-01  1.2e-02  1.3e-02  6.79e-01   4.904313668e+02   4.905719210e+02   1.2e-02  97.57 
    8   3.7e-01  9.1e-03  1.0e-02  6.18e-01   4.944850631e+02   4.946792913e+02   9.1e-03  98.12 
    9   2.5e-01  6.1e-03  7.2e-03  2.96e-01   6.129255413e+02   6.131604659e+02   6.1e-03  98.68 
    10  1.6e-01  3.8e-03  4.5e-03  2.50e-01   8.308572820e+02   8.311056790e+02   3.8e-03  99.22 
    11  9.8e-02  2.4e-03  2.3e-03  6.55e-01   1.139192306e+03   1.139369915e+03   2.4e-03  99.77 
    12  7.2e-02  1.8e-03  1.7e-03  4.11e-01   1.303282105e+03   1.303457673e+03   1.8e-03  100.32
    13  2.2e-02  5.3e-04  3.8e-04  4.83e-01   1.844543644e+03   1.844648722e+03   5.3e-04  100.89
    14  6.5e-03  1.6e-04  6.6e-05  8.74e-01   2.147120530e+03   2.147157160e+03   1.6e-04  101.51
    15  6.7e-04  1.6e-05  2.6e-06  8.99e-01   2.291660423e+03   2.291666225e+03   1.6e-05  102.29
    16  1.3e-04  3.2e-06  2.3e-07  9.94e-01   2.309537827e+03   2.309539007e+03   3.2e-06  102.87
    17  4.6e-05  1.1e-06  4.8e-08  9.88e-01   2.312484238e+03   2.312484647e+03   1.1e-06  103.54
    18  1.8e-05  4.3e-07  1.2e-08  9.98e-01   2.313457331e+03   2.313457490e+03   4.3e-07  104.74
    19  9.5e-07  2.3e-08  6.8e-11  9.99e-01   2.314036746e+03   2.314036756e+03   2.3e-08  105.74
    20  9.5e-07  2.3e-08  6.8e-11  1.00e+00   2.314036746e+03   2.314036756e+03   2.3e-08  106.96
    21  9.5e-07  2.3e-08  6.8e-11  1.00e+00   2.314036746e+03   2.314036756e+03   2.3e-08  108.25
    Optimizer terminated. Time: 110.19  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 2.3140367463e+03    nrm: 5e+03    Viol.  con: 5e-07    var: 0e+00    cones: 0e+00  
      Dual.    obj: 2.3140367563e+03    nrm: 1e+05    Viol.  con: 1e-03    var: 3e-08    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)








.. image:: /auto_examples/coil_design/images/sphx_glr_biplanar_coil_design_001.png
    :class: sphx-glr-single-img




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





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 5 minutes  43.006 seconds)

**Estimated memory usage:**  4591 MB


.. _sphx_glr_download_auto_examples_coil_design_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: biplanar_coil_design.py <biplanar_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: biplanar_coil_design.ipynb <biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
