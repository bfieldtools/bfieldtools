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

    from bfieldtools.mesh_class import Conductor
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

    #Create Conductor object, which finds the holes and sets the boundary condition
    coil = Conductor(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)







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

    Computing magnetic field coupling matrix, 2772 vertices by 160 target points... took 0.23 seconds.



Run QP solver


.. code-block:: default

    import mosek

    coil.s, prob = optimize_streamfunctions(coil,
                                       bfield_specification,
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 27549 MiB required for 2772 by 2772 vertices...
    Computing inductance matrix in 60 chunks (10101 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 38.12 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 3364            
      Cones                  : 1               
      Scalar variables       : 5177            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 3364            
      Cones                  : 1               
      Scalar variables       : 5177            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2403
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 3364              conic                  : 2404            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.73              dense det. time        : 0.00            
    Factor     - ML order time          : 0.17              GP order time          : 0.00            
    Factor     - nonzeros before factor : 2.89e+06          after factor           : 2.89e+06        
    Factor     - dense dim.             : 0                 flops                  : 2.13e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  48.81 
    1   3.6e+01  5.6e-01  1.5e+00  -6.20e-01  1.057039411e+01   1.024347864e+01   5.6e-01  49.08 
    2   3.2e+01  5.1e-01  2.3e-01  8.02e-01   3.529324593e+01   3.427460253e+01   5.1e-01  49.33 
    3   1.4e+01  2.1e-01  1.4e-01  -1.12e-01  7.399485968e+01   7.350916188e+01   2.1e-01  49.59 
    4   6.7e+00  1.0e-01  1.3e-01  2.45e-01   7.816471843e+01   7.777292693e+01   1.0e-01  49.84 
    5   2.1e+00  3.3e-02  2.9e-02  4.36e-01   1.951267296e+02   1.950674157e+02   3.3e-02  50.10 
    6   7.0e-01  1.1e-02  7.2e-03  7.37e-01   3.333579805e+02   3.333701148e+02   1.1e-02  50.39 
    7   3.0e-01  4.7e-03  2.3e-03  9.23e-01   4.014258016e+02   4.014427924e+02   4.7e-03  50.64 
    8   3.2e-02  5.1e-04  1.0e-04  9.58e-01   4.507832760e+02   4.507887089e+02   5.1e-04  50.93 
    9   4.2e-03  6.5e-05  4.7e-06  1.03e+00   4.570326404e+02   4.570333241e+02   6.5e-05  51.22 
    10  1.1e-03  1.7e-05  6.3e-07  1.01e+00   4.577867972e+02   4.577869831e+02   1.7e-05  51.47 
    11  9.6e-05  1.5e-06  1.6e-08  1.00e+00   4.580310515e+02   4.580310678e+02   1.5e-06  51.77 
    12  1.5e-05  2.4e-07  1.0e-09  1.00e+00   4.580507069e+02   4.580507095e+02   2.4e-07  52.04 
    13  3.0e-07  6.5e-08  3.7e-12  1.00e+00   4.580543203e+02   4.580543203e+02   4.7e-09  52.39 
    14  1.6e-07  1.3e-08  1.1e-12  1.00e+00   4.580543787e+02   4.580543787e+02   9.3e-10  52.90 
    15  6.6e-07  1.5e-09  9.9e-14  1.00e+00   4.580543915e+02   4.580543915e+02   1.1e-10  53.33 
    Optimizer terminated. Time: 53.53   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 4.5805439148e+02    nrm: 9e+02    Viol.  con: 5e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 4.5805439149e+02    nrm: 1e+04    Viol.  con: 0e+00    var: 3e-09    cones: 0e+00  



Plot the computed streamfunction


.. code-block:: default


    coil.s.plot(ncolors=256)



.. image:: /auto_examples/coil_design/images/sphx_glr_coil_with_holes_001.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  13.783 seconds)

**Estimated memory usage:**  2662 MB


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
