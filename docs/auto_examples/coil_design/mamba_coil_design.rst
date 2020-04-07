.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_coil_design_mamba_coil_design.py>`     to download the full example code
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
            "bfieldtools", "example_meshes/10x10_plane_hires.obj"
        ),
        process=False,
    )

    planemesh.apply_scale(scaling_factor)

    # planemesh.vertices, planemesh.faces = trimesh.remesh.subdivide(planemesh.vertices, planemesh.faces)


    # Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 1.5, 0]) * scaling_factor

    # Create coil plane pairs
    coil_plus = trimesh.Trimesh(
        planemesh.vertices + center_offset + standoff, planemesh.faces, process=False
    )

    coil_minus = trimesh.Trimesh(
        planemesh.vertices + center_offset - standoff, planemesh.faces, process=False
    )

    joined_planes = coil_plus.union(coil_minus)

    # Create mesh class object
    coil = Conductor(
        verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True
    )








Set up target and stray field points. Here, the target points are on a planar
4x4 grid slightly smaller than the coil dimensions.


.. code-block:: default


    center = np.array([0, 0, 0]) * scaling_factor

    sidelength = 0.5 * scaling_factor
    n = 4

    height = 0.1
    n_height = 2
    xx = np.linspace(-sidelength / 2, sidelength / 2, n)
    yy = np.linspace(-height / 2, height / 2, n_height)
    zz = np.linspace(-sidelength / 2, sidelength / 2, n)
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing="ij")

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

    target_points = np.asarray(grid_target_points).reshape((-1, 3))
    target_field = np.asarray(target_field).reshape((-1,))

    target_field = np.array(
        [np.zeros((len(target_field),)), target_field, np.zeros((len(target_field),))]
    ).T


    target_rel_error = np.zeros_like(target_field)
    target_rel_error[:, 1] += 0.05

    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 1] += 0.01
    target_abs_error[:, 0::2] += 0.05








Plot target points and mesh


.. code-block:: default

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))

    mlab.quiver3d(*target_points.T, *target_field.T)
    coil.plot_mesh()





.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/coil_design/images/sphx_glr_mamba_coil_design_001.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/coil_design/images/sphx_glr_mamba_coil_design_002.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.surface.Surface object at 0x00000254536BC830>



Compute coupling matrix that is used to compute the generated magnetic field, create field specification


.. code-block:: default



    target_spec = {
        "coupling": coil.B_coupling(target_points),
        "rel_error": target_rel_error,
        "abs_error": target_abs_error,
        "target": target_field,
    }





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 512 target points... took 1.09 seconds.




Run QP solver


.. code-block:: default


    import mosek

    coil.j, prob = optimize_streamfunctions(
        coil,
        [target_spec],
        objective="minimum_inductive_energy",
        solver="MOSEK",
        solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
    )






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 100 chunks (8384 MiB memory free),                  when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 40.34 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


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
    Factor     - setup time             : 1.38              dense det. time        : 0.00            
    Factor     - ML order time          : 0.27              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.53e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   2.5e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  132.41
    1   1.1e+01  4.5e-01  4.3e-01  3.73e-01   9.145114752e+01   9.082987855e+01   4.5e-01  133.92
    2   3.4e+00  1.3e-01  7.2e-02  7.00e-01   1.635630933e+02   1.633665156e+02   1.3e-01  135.39
    3   1.4e+00  5.7e-02  2.5e-02  9.61e-01   1.755024270e+02   1.754262626e+02   5.7e-02  136.84
    4   2.8e-02  1.1e-03  4.9e-05  1.13e+00   1.893939620e+02   1.893923008e+02   1.1e-03  138.52
    5   7.5e-03  2.9e-04  8.4e-06  1.01e+00   1.894767240e+02   1.894763509e+02   2.9e-04  140.00
    6   5.8e-04  2.3e-05  2.1e-07  1.00e+00   1.895531186e+02   1.895530951e+02   2.3e-05  141.42
    7   5.3e-05  2.1e-06  5.9e-09  1.00e+00   1.895628473e+02   1.895628452e+02   2.1e-06  142.97
    8   2.7e-05  1.1e-06  2.1e-09  1.00e+00   1.895633339e+02   1.895633328e+02   1.1e-06  144.39
    9   1.0e-05  4.0e-07  5.2e-10  1.00e+00   1.895636573e+02   1.895636569e+02   4.0e-07  145.92
    10  1.4e-06  5.5e-08  2.4e-11  1.00e+00   1.895638193e+02   1.895638192e+02   5.5e-08  147.48
    11  2.0e-07  2.1e-09  1.6e-12  1.00e+00   1.895638461e+02   1.895638462e+02   1.1e-10  149.11
    Optimizer terminated. Time: 149.66  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.8956384607e+02    nrm: 4e+02    Viol.  con: 3e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.8956384623e+02    nrm: 1e+02    Viol.  con: 5e-10    var: 5e-10    cones: 0e+00  




Plot coil windings and target points


.. code-block:: default


    loops = scalar_contour(coil.mesh, coil.j, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors="auto", figure=f, tube_radius=0.025)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)

    f.scene.isometric_view()



.. image:: /auto_examples/coil_design/images/sphx_glr_mamba_coil_design_003.png
    :class: sphx-glr-single-img






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 4 minutes  16.597 seconds)


.. _sphx_glr_download_auto_examples_coil_design_mamba_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: mamba_coil_design.py <mamba_coil_design.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: mamba_coil_design.ipynb <mamba_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
