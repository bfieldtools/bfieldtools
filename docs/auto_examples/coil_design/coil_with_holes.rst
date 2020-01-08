.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_coil_with_holes.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_coil_with_holes.py:


Coil with interior holes
========================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes. The coil planes have holes in them,



.. code-block:: default


    import numpy as np
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
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/plane_w_holes.stl'), process=False)

    angle=np.pi/2
    rotation_matrix = np.array([[1, 0, 0, 0],
                                [0, np.cos(angle), -np.sin(angle), 0],
                                [0, np.sin(angle), np.cos(angle), 0],
                                [0, 0, 0, 1]
                                  ])

    planemesh.apply_transform(rotation_matrix)
    planemesh.apply_scale(scaling_factor)

    #Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 20, 0]) * scaling_factor

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

    sidelength = 10 * scaling_factor
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










Let's find and separate the inner and outer boundaries of the coil mesh


.. code-block:: default


    inner_bounds = np.intersect1d(coil.boundary_verts, np.where(np.linalg.norm(coil.mesh.vertices[:,0::2], axis=1)< 0.015)[0])

    centre_hole1 = np.intersect1d(np.intersect1d(coil.boundary_verts,
                                                 np.where(np.linalg.norm(coil.mesh.vertices[:,0::2], axis=1)< 0.004)[0]),
                                  np.where(coil.mesh.vertices[:,1] < 0)[0])

    centre_hole2 = np.intersect1d(np.intersect1d(coil.boundary_verts,
                                                 np.where(np.linalg.norm(coil.mesh.vertices[:,0::2], axis=1)< 0.004)[0]),
                                  np.where(coil.mesh.vertices[:,1] > 0)[0])

    left_hole1 = np.intersect1d(np.intersect1d(np.intersect1d(coil.boundary_verts,
                                               np.where(coil.mesh.vertices[:,0] < -0.004)[0]), inner_bounds),
                                np.where(coil.mesh.vertices[:,1] < 0)[0])

    left_hole2 = np.intersect1d(np.intersect1d(np.intersect1d(coil.boundary_verts,
                                               np.where(coil.mesh.vertices[:,0] < -0.004)[0]), inner_bounds),
                                np.where(coil.mesh.vertices[:,1] > 0)[0])

    right_hole1 = np.intersect1d(np.intersect1d(np.intersect1d(coil.boundary_verts,
                                               np.where(coil.mesh.vertices[:,0] > 0.004)[0]), inner_bounds),
                                np.where(coil.mesh.vertices[:,1] < 0)[0])

    right_hole2 = np.intersect1d(np.intersect1d(np.intersect1d(coil.boundary_verts,
                                               np.where(coil.mesh.vertices[:,0] > 0.004)[0]), inner_bounds),
                                np.where(coil.mesh.vertices[:,1] > 0)[0])

    outer_bounds = np.setdiff1d(coil.boundary_verts, inner_bounds)


    graph = trimesh.graph.vertex_adjacency_graph(coil.mesh)

    zero_eq_indices = outer_bounds
    iso_eq_indices = [left_hole1, centre_hole1, right_hole1, left_hole2, centre_hole2, right_hole2]

    boundary_constraints = {'zero_eq_indices':zero_eq_indices , 'iso_eq_indices': iso_eq_indices}







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

    bfield_specification = [target_spec]





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 2772 vertices by 160 target points... took 0.25 seconds.



Run QP solver


.. code-block:: default

    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       bfield_specification,
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}},
                                       boundary_constraints=boundary_constraints
                                       )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 307359 MiB required for 2772 times 2772 vertices...
    Computing inductance matrix in 32 chunks since 9832 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 60.82 seconds.
    Pre-existing problem not passed, creating...
    Passing boundary zero equality constraint
    Passing boundary equality constraints
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 4108            
      Cones                  : 1               
      Scalar variables       : 5545            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 4108            
      Cones                  : 1               
      Scalar variables       : 5545            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2397
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 3732              conic                  : 2772            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.70              dense det. time        : 0.00            
    Factor     - ML order time          : 0.10              GP order time          : 0.00            
    Factor     - nonzeros before factor : 2.87e+06          after factor           : 2.87e+06        
    Factor     - dense dim.             : 0                 flops                  : 2.33e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  52.54 
    1   6.0e+01  4.7e-01  1.3e+00  -8.27e-01  1.384962869e+01   1.386946830e+01   4.7e-01  52.84 
    2   5.1e+01  4.0e-01  8.6e-01  -2.18e-01  7.010307566e+01   6.976387568e+01   4.0e-01  53.10 
    3   2.1e+01  1.7e-01  5.5e-01  -6.24e-01  1.292869617e+02   1.304411776e+02   1.7e-01  53.45 
    4   1.3e+01  1.0e-01  1.6e-01  6.79e-01   4.226105729e+02   4.224895934e+02   1.0e-01  53.72 
    5   4.3e+00  3.4e-02  8.4e-02  2.63e-03   6.636826144e+02   6.645554835e+02   3.4e-02  54.00 
    6   8.0e-01  6.2e-03  6.4e-03  1.63e-01   1.668116368e+03   1.668209765e+03   6.2e-03  54.37 
    7   3.7e-01  2.8e-03  2.9e-03  6.84e-01   1.984140844e+03   1.984310583e+03   2.8e-03  54.64 
    8   4.8e-02  3.7e-04  2.2e-04  9.02e-01   2.299257269e+03   2.299333600e+03   3.7e-04  54.98 
    9   4.3e-03  3.3e-05  6.0e-06  1.03e+00   2.356815098e+03   2.356822154e+03   3.3e-05  55.33 
    10  1.3e-03  1.0e-05  1.1e-06  1.00e+00   2.361140207e+03   2.361142487e+03   1.0e-05  55.61 
    11  1.3e-04  1.0e-06  3.2e-08  9.97e-01   2.362979138e+03   2.362979364e+03   1.0e-06  55.91 
    12  1.5e-07  1.1e-07  8.8e-12  1.00e+00   2.363174627e+03   2.363174628e+03   1.6e-09  56.27 
    13  2.0e-06  5.7e-08  1.3e-11  1.00e+00   2.363174737e+03   2.363174738e+03   8.0e-10  56.72 
    14  2.2e-06  4.3e-08  1.4e-11  1.00e+00   2.363174765e+03   2.363174765e+03   6.0e-10  57.20 
    15  2.2e-06  4.3e-08  1.4e-11  1.00e+00   2.363174765e+03   2.363174765e+03   6.0e-10  57.75 
    16  2.2e-06  4.3e-08  1.4e-11  1.00e+00   2.363174765e+03   2.363174765e+03   6.0e-10  58.33 
    Optimizer terminated. Time: 59.13   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 2.3631747647e+03    nrm: 5e+03    Viol.  con: 7e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 2.3631747655e+03    nrm: 6e+04    Viol.  con: 0e+00    var: 2e-07    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    N_contours = 20

    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.1)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)


.. image:: /auto_examples/coil_design/images/sphx_glr_coil_with_holes_001.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  39.009 seconds)

**Estimated memory usage:**  2746 MB


.. _sphx_glr_download_auto_examples_coil_design_coil_with_holes.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: coil_with_holes.py <coil_with_holes.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: coil_with_holes.ipynb <coil_with_holes.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
