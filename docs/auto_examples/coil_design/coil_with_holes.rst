.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_coil_design_coil_with_holes.py>`     to download the full example code
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

    from bfieldtools.conductor import Conductor
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

    Computing magnetic field coupling matrix, 2772 vertices by 160 target points... took 0.29 seconds.




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
    Computing inductance matrix in 80 chunks (8096 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 51.41 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 3732            
      Cones                  : 1               
      Scalar variables       : 5545            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 3732            
      Cones                  : 1               
      Scalar variables       : 5545            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2773
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 3732              conic                  : 2772            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.75              dense det. time        : 0.00            
    Factor     - ML order time          : 0.17              GP order time          : 0.00            
    Factor     - nonzeros before factor : 3.85e+06          after factor           : 3.85e+06        
    Factor     - dense dim.             : 0                 flops                  : 3.21e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   3.2e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  82.08 
    1   2.2e+01  6.9e-01  9.8e-01  3.92e-02   3.841820107e+00   3.041284579e+00   6.9e-01  83.88 
    2   1.5e+01  4.8e-01  1.2e-01  2.65e-01   1.450577318e+01   1.394852532e+01   4.8e-01  84.78 
    3   7.6e+00  2.4e-01  5.5e-02  2.40e+00   2.470331119e+01   2.455714325e+01   2.4e-01  85.70 
    4   1.3e+00  4.0e-02  5.4e-03  1.73e+00   2.666954569e+01   2.665266899e+01   4.0e-02  86.69 
    5   3.4e-01  1.1e-02  7.1e-04  1.08e+00   2.728123817e+01   2.727662870e+01   1.1e-02  87.72 
    6   2.5e-01  7.9e-03  5.3e-04  6.43e-01   2.732205442e+01   2.731806651e+01   7.9e-03  88.59 
    7   1.4e-01  4.5e-03  2.2e-04  7.70e-01   2.746821574e+01   2.746581829e+01   4.5e-03  89.47 
    8   2.5e-02  7.8e-04  1.4e-05  5.94e-01   2.766833285e+01   2.766780009e+01   7.8e-04  90.39 
    9   8.7e-05  2.7e-06  3.0e-09  9.43e-01   2.771597599e+01   2.771597409e+01   2.7e-06  91.28 
    10  8.3e-07  3.5e-08  3.1e-12  1.00e+00   2.771616534e+01   2.771616532e+01   2.6e-08  93.97 
    Optimizer terminated. Time: 94.23   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 2.7716165339e+01    nrm: 6e+01    Viol.  con: 1e-08    var: 0e+00    cones: 0e+00  
      Dual.    obj: 2.7716165318e+01    nrm: 3e+02    Viol.  con: 2e-08    var: 2e-09    cones: 0e+00  




Plot the computed streamfunction


.. code-block:: default


    coil.s.plot(ncolors=256)



.. image:: /auto_examples/coil_design/images/sphx_glr_coil_with_holes_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.surface.Surface object at 0x000001504D6BE9E8>




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  51.801 seconds)


.. _sphx_glr_download_auto_examples_coil_design_coil_with_holes.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: coil_with_holes.py <coil_with_holes.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: coil_with_holes.ipynb <coil_with_holes.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
