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

    Computing C matrix, 2772 vertices by 160 target points... took 0.28 seconds.



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

    Computing inductance matrix in 2 chunks since 7 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 52.12 seconds.
    Pre-existing problem not passed, creating...
    Passing boundary zero equality constraint
    Passing boundary equality constraints
    Passing parameters to problem...
    Passing problem to solver...


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
    Factor     - setup time             : 0.78              dense det. time        : 0.00            
    Factor     - ML order time          : 0.17              GP order time          : 0.00            
    Factor     - nonzeros before factor : 2.87e+06          after factor           : 2.87e+06        
    Factor     - dense dim.             : 0                 flops                  : 2.33e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   2.6e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  47.18 
    1   1.3e+02  4.9e-01  1.2e+00  -8.01e-01  1.314016891e+01   1.297374658e+01   4.9e-01  47.47 
    2   7.2e+01  2.8e-01  8.5e-01  -7.16e-01  5.227631276e+01   5.293911450e+01   2.8e-01  47.74 
    3   4.6e+01  1.8e-01  3.6e-01  2.31e-01   2.416119265e+02   2.416233728e+02   1.8e-01  48.02 
    4   1.7e+01  6.7e-02  2.1e-01  -3.74e-01  3.487844131e+02   3.502915686e+02   6.7e-02  48.29 
    5   6.0e+00  2.3e-02  4.5e-02  2.30e-01   1.041175646e+03   1.041708335e+03   2.3e-02  48.56 
    6   2.1e+00  8.1e-03  1.6e-02  1.88e-01   1.654937539e+03   1.655654721e+03   8.1e-03  48.84 
    7   4.8e-01  1.9e-03  1.7e-03  9.81e-01   2.165332079e+03   2.165485403e+03   1.9e-03  49.12 
    8   4.3e-03  1.7e-05  9.7e-07  9.78e-01   2.331560695e+03   2.331561099e+03   1.7e-05  49.49 
    9   1.6e-05  6.1e-08  1.6e-10  1.00e+00   2.333080133e+03   2.333080133e+03   6.1e-08  49.87 
    10  1.6e-07  3.5e-10  3.2e-12  1.00e+00   2.333086102e+03   2.333086103e+03   7.9e-11  50.47 
    Optimizer terminated. Time: 50.69   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 2.3330861020e+03    nrm: 5e+03    Viol.  con: 1e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 2.3330861028e+03    nrm: 6e+04    Viol.  con: 1e-07    var: 1e-09    cones: 0e+00  



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

   **Total running time of the script:** ( 2 minutes  22.594 seconds)

**Estimated memory usage:**  5621 MB


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
