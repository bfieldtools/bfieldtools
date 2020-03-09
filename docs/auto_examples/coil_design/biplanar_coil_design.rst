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

    from bfieldtools.mesh_class import Conductor
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
    coil = Conductor(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)







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

    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.28 seconds.
    Computing magnetic field coupling matrix, 3184 vertices by 642 target points... took 0.70 seconds.



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

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 60 chunks (11699 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 43.59 seconds.
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
      Scalar variables       : 6083            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 7710            
      Cones                  : 1               
      Scalar variables       : 6083            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 7710              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 1.88              dense det. time        : 0.00            
    Factor     - ML order time          : 0.17              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 5.26e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  90.59 
    1   2.8e+01  4.4e-01  1.2e+00  -8.72e-01  7.505693341e+00   7.568296918e+00   4.4e-01  91.24 
    2   1.9e+01  2.9e-01  8.9e-01  -6.59e-01  3.189700433e+01   3.251077617e+01   2.9e-01  91.79 
    3   1.1e+01  1.7e-01  5.4e-01  -5.00e-01  1.413946429e+02   1.426196877e+02   1.7e-01  92.34 
    4   7.2e+00  1.1e-01  3.6e-01  -2.45e-01  1.782403672e+02   1.796636046e+02   1.1e-01  92.89 
    5   2.7e+00  4.2e-02  1.2e-01  -1.72e-02  4.129910060e+02   4.142915941e+02   4.2e-02  93.44 
    6   9.1e-01  1.4e-02  2.6e-02  4.89e-01   5.064346055e+02   5.070514328e+02   1.4e-02  94.11 
    7   7.5e-01  1.2e-02  2.2e-02  5.44e-01   5.373466093e+02   5.380158842e+02   1.2e-02  94.66 
    8   5.6e-01  8.8e-03  1.6e-02  4.81e-01   6.793322975e+02   6.799762595e+02   8.8e-03  95.22 
    9   3.5e-01  5.4e-03  8.1e-03  6.99e-01   1.057425454e+03   1.057850456e+03   5.4e-03  95.82 
    10  2.3e-01  3.7e-03  5.4e-03  4.88e-01   1.177304957e+03   1.177741132e+03   3.7e-03  96.39 
    11  5.8e-02  9.0e-04  7.7e-04  8.60e-01   1.853290108e+03   1.853444126e+03   9.0e-04  97.13 
    12  1.3e-02  2.1e-04  9.3e-05  8.99e-01   2.192215741e+03   2.192260123e+03   2.1e-04  97.78 
    13  3.7e-04  5.7e-06  4.6e-07  9.58e-01   2.310530307e+03   2.310531694e+03   5.7e-06  98.56 
    14  2.1e-06  3.3e-08  8.2e-11  9.99e-01   2.314263581e+03   2.314263590e+03   3.3e-08  99.32 
    15  2.0e-06  3.1e-08  8.2e-11  1.00e+00   2.314264727e+03   2.314264736e+03   3.1e-08  100.46
    16  2.0e-06  3.1e-08  9.1e-11  9.99e-01   2.314264738e+03   2.314264747e+03   3.1e-08  101.60
    17  9.4e-07  1.5e-08  1.5e-11  1.00e+00   2.314275954e+03   2.314275958e+03   1.5e-08  102.66
    18  7.0e-07  1.1e-08  5.7e-11  1.00e+00   2.314278463e+03   2.314278466e+03   1.1e-08  103.73
    19  4.0e-07  6.2e-09  1.4e-11  1.00e+00   2.314281761e+03   2.314281763e+03   6.2e-09  104.80
    20  2.1e-07  3.2e-09  2.9e-12  1.00e+00   2.314283819e+03   2.314283820e+03   3.2e-09  105.87
    21  1.0e-07  1.6e-09  1.5e-12  1.00e+00   2.314284931e+03   2.314284931e+03   1.6e-09  106.95
    22  2.1e-07  8.2e-10  5.0e-12  1.00e+00   2.314285499e+03   2.314285499e+03   8.2e-10  108.03
    23  3.3e-08  4.1e-10  1.5e-12  1.00e+00   2.314285783e+03   2.314285783e+03   4.1e-10  109.10
    Optimizer terminated. Time: 109.63  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 2.3142857834e+03    nrm: 5e+03    Viol.  con: 7e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 2.3142857835e+03    nrm: 1e+05    Viol.  con: 2e-05    var: 2e-09    cones: 0e+00  



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

   **Total running time of the script:** ( 5 minutes  39.392 seconds)

**Estimated memory usage:**  4685 MB


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
