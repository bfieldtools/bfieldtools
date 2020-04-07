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


    # Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    # Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(
        file_obj=pkg_resources.resource_filename(
            "bfieldtools", "example_meshes/plane_w_holes.stl"
        ),
        process=False,
    )

    angle = np.pi / 2
    rotation_matrix = np.array(
        [
            [1, 0, 0, 0],
            [0, np.cos(angle), -np.sin(angle), 0],
            [0, np.sin(angle), np.cos(angle), 0],
            [0, 0, 0, 1],
        ]
    )

    planemesh.apply_transform(rotation_matrix)
    planemesh.apply_scale(scaling_factor)

    # Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 20, 0]) * scaling_factor

    # Create coil plane pairs
    coil_plus = trimesh.Trimesh(
        planemesh.vertices + center_offset + standoff, planemesh.faces, process=False
    )

    coil_minus = trimesh.Trimesh(
        planemesh.vertices + center_offset - standoff, planemesh.faces, process=False
    )

    joined_planes = coil_plus.union(coil_minus)

    # Create Conductor object, which finds the holes and sets the boundary condition
    coil = Conductor(
        verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True
    )








Set up target and stray field points


.. code-block:: default


    # Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 0]) * scaling_factor

    sidelength = 10 * scaling_factor
    n = 8
    xx = np.linspace(-sidelength / 2, sidelength / 2, n)
    yy = np.linspace(-sidelength / 2, sidelength / 2, n)
    zz = np.linspace(-sidelength / 2, sidelength / 2, n)
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing="ij")

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    target_points = np.array([x, y, z]).T

    # Turn cube into sphere by rejecting points "in the corners"
    target_points = (
        target_points[np.linalg.norm(target_points, axis=1) < sidelength / 2] + center
    )









Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    # The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function

    target_field = np.zeros(target_points.shape)
    target_field[:, 0] = target_field[:, 0] + 1

    target_rel_error = np.zeros_like(target_field)
    target_rel_error[:, 0] += 0.01

    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 0] += 0.001
    target_abs_error[:, 1:3] += 0.005

    target_spec = {
        "coupling": coil.B_coupling(target_points),
        "rel_error": target_rel_error,
        "abs_error": target_abs_error,
        "target": target_field,
    }

    bfield_specification = [target_spec]





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 2772 vertices by 160 target points... took 0.33 seconds.




Run QP solver


.. code-block:: default

    import mosek

    coil.s, prob = optimize_streamfunctions(
        coil,
        bfield_specification,
        objective="minimum_inductive_energy",
        solver="MOSEK",
        solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
    )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Estimating 27549 MiB required for 2772 by 2772 vertices...
    Computing inductance matrix in 80 chunks (8998 MiB memory free),                  when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 25.46 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 3364            
      Cones                  : 1               
      Scalar variables       : 4807            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 3364            
      Cones                  : 1               
      Scalar variables       : 4807            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2403
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 3364              conic                  : 2404            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.63              dense det. time        : 0.00            
    Factor     - ML order time          : 0.16              GP order time          : 0.00            
    Factor     - nonzeros before factor : 2.89e+06          after factor           : 2.89e+06        
    Factor     - dense dim.             : 0                 flops                  : 2.13e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   3.2e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  67.17 
    1   2.2e+01  6.8e-01  6.0e-01  5.66e-01   8.753280039e+00   7.935396740e+00   6.8e-01  67.89 
    2   1.3e+01  4.2e-01  2.3e-01  7.70e-01   2.207612917e+01   2.163215677e+01   4.2e-01  68.59 
    3   5.6e+00  1.7e-01  9.2e-02  2.08e+00   2.871617786e+01   2.857708476e+01   1.7e-01  69.39 
    4   2.6e+00  8.0e-02  3.3e-02  7.19e-01   3.872879099e+01   3.865848781e+01   8.0e-02  70.14 
    5   5.3e-01  1.6e-02  2.4e-03  9.32e-01   4.848342223e+01   4.846630601e+01   1.6e-02  70.83 
    6   9.0e-02  2.8e-03  1.7e-04  9.85e-01   5.037994168e+01   5.037689232e+01   2.8e-03  71.59 
    7   6.6e-03  2.1e-04  4.1e-06  9.06e-01   5.088983144e+01   5.088962150e+01   2.1e-04  72.28 
    8   2.7e-04  8.4e-06  3.6e-08  9.92e-01   5.094421871e+01   5.094421067e+01   8.4e-06  72.97 
    9   5.9e-05  2.7e-06  4.1e-09  1.00e+00   5.094641132e+01   5.094640967e+01   2.0e-06  73.67 
    10  4.6e-06  2.1e-07  1.0e-10  1.00e+00   5.094699440e+01   5.094699423e+01   1.6e-07  74.38 
    11  1.0e-07  1.9e-09  1.0e-13  1.00e+00   5.094704316e+01   5.094704317e+01   1.0e-09  75.11 
    Optimizer terminated. Time: 75.34   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 5.0947043165e+01    nrm: 1e+02    Viol.  con: 7e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 5.0947043168e+01    nrm: 6e+02    Viol.  con: 0e+00    var: 7e-10    cones: 0e+00  




Plot the computed streamfunction


.. code-block:: default


    coil.s.plot(ncolors=256)



.. image:: /auto_examples/coil_design/images/sphx_glr_coil_with_holes_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.surface.Surface object at 0x0000025453693468>




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  4.154 seconds)


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
