.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_mamba_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_mamba_coil_design.py:


MAMBA coil
==========

Compact example of a biplanar coil producing homogeneous field in a number of target
regions arranged in a grid.



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
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

    planemesh.apply_scale(scaling_factor)

    #planemesh.vertices, planemesh.faces = trimesh.remesh.subdivide(planemesh.vertices, planemesh.faces)


    #Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 1.5, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                             planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                         planemesh.faces, process=False)

    joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = MeshWrapper(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)







Set up target and stray field points. Here, the target points are on a planar
4x4 grid slightly smaller than the coil dimensions.


.. code-block:: default


    center = np.array([0, 0, 0]) * scaling_factor

    sidelength = 0.5 * scaling_factor
    n = 4

    height = 0.1
    n_height = 2
    xx = np.linspace(-sidelength/2, sidelength/2, n)
    yy = np.linspace(-height/2, height/2, n_height)
    zz = np.linspace(-sidelength/2, sidelength/2, n)
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    target_points = np.array([x, y, z]).T


    grid_target_points = list()
    target_field = list()

    hori_offsets = [-3, -1, 1, 3]
    vert_offsets = [-3, -1, 1, 3]

    for i, offset_x in enumerate(hori_offsets):
        for j, offset_y in enumerate(vert_offsets):
            grid_target_points.append(target_points + np.array([offset_x, 0, offset_y]))
            target_field.append((i + j - 3) * np.ones((len(target_points),)))

    target_points = np.asarray(grid_target_points).reshape((-1,3))
    target_field = np.asarray(target_field).reshape((-1,))

    target_field = np.array([np.zeros((len(target_field),)), target_field, np.zeros((len(target_field),))]).T


    target_rel_error = np.zeros_like(target_field)
    target_rel_error[:, 1] += 0.05

    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 1] += 0.01
    target_abs_error[:, 0::2] += 0.05







Plot target points and mesh


.. code-block:: default

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(800, 800))

    mlab.quiver3d(*target_points.T, *target_field.T)
    coil.plot_mesh()





.. image:: /auto_examples/coil_design/images/sphx_glr_mamba_coil_design_001.png
    :class: sphx-glr-single-img




Compute coupling matrix that is used to compute the generated magnetic field, create field specification


.. code-block:: default



    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 512 target points... took 0.68 seconds.



Run QP solver


.. code-block:: default


    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 405514 MiB required for 3184 times 3184 vertices...
    Computing inductance matrix in 42 chunks since 9800 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 86.23 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 5970            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 5970            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 5970              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 1.62              dense det. time        : 0.00            
    Factor     - ML order time          : 0.30              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.53e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   2.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  92.49 
    1   1.1e+01  4.6e-01  4.8e-01  2.12e-01   6.597411176e+01   6.533623471e+01   4.6e-01  93.07 
    2   6.0e+00  2.5e-01  2.3e-01  2.10e-01   1.540986144e+02   1.537101559e+02   2.5e-01  93.67 
    3   5.0e+00  2.1e-01  1.9e-01  6.31e-01   1.799383624e+02   1.796052015e+02   2.1e-01  94.17 
    4   9.7e-01  4.0e-02  1.7e-02  7.68e-01   3.380757286e+02   3.380076296e+02   4.0e-02  94.86 
    5   5.0e-01  2.1e-02  7.3e-03  8.87e-01   3.735298693e+02   3.734982875e+02   2.1e-02  95.40 
    6   2.9e-01  1.2e-02  3.5e-03  9.19e-01   3.945710275e+02   3.945537579e+02   1.2e-02  95.91 
    7   1.4e-01  5.7e-03  1.2e-03  9.36e-01   4.128020623e+02   4.127949103e+02   5.7e-03  96.48 
    8   2.1e-02  8.8e-04  8.3e-05  9.65e-01   4.284418731e+02   4.284412976e+02   8.8e-04  97.07 
    9   2.7e-03  1.1e-04  3.8e-06  9.95e-01   4.319789642e+02   4.319788943e+02   1.1e-04  97.60 
    10  2.6e-05  1.1e-06  3.7e-09  9.97e-01   4.325195035e+02   4.325195028e+02   1.1e-06  98.17 
    11  3.4e-06  1.4e-07  2.0e-10  1.00e+00   4.325241783e+02   4.325241781e+02   1.4e-07  98.73 
    12  2.7e-07  6.6e-09  4.5e-12  1.00e+00   4.325248512e+02   4.325248513e+02   6.6e-09  99.85 
    Optimizer terminated. Time: 100.27  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 4.3252485119e+02    nrm: 9e+02    Viol.  con: 3e-08    var: 0e+00    cones: 0e+00  
      Dual.    obj: 4.3252485126e+02    nrm: 5e+02    Viol.  con: 0e+00    var: 2e-09    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.025)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)

    f.scene.isometric_view()



.. image:: /auto_examples/coil_design/images/sphx_glr_mamba_coil_design_002.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 4 minutes  13.713 seconds)

**Estimated memory usage:**  3678 MB


.. _sphx_glr_download_auto_examples_coil_design_mamba_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: mamba_coil_design.py <mamba_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: mamba_coil_design.ipynb <mamba_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
