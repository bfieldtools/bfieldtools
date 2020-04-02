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
    coil = Conductor(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)








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


    <mayavi.modules.surface.Surface object at 0x000001B5C4CB2258>



Compute coupling matrix that is used to compute the generated magnetic field, create field specification


.. code-block:: default



    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 512 target points... took 0.92 seconds.




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

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 160 chunks (4872 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 59.93 seconds.
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
    Factor     - setup time             : 0.97              dense det. time        : 0.00            
    Factor     - ML order time          : 0.17              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.53e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   5.1e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  110.84
    1   2.3e+01  4.5e-01  4.3e-01  3.68e-01   9.097536404e+01   9.035380857e+01   4.5e-01  112.11
    2   6.7e+00  1.3e-01  7.2e-02  6.98e-01   1.640403677e+02   1.638441069e+02   1.3e-01  113.41
    3   2.9e+00  5.7e-02  2.5e-02  9.63e-01   1.756176942e+02   1.755420942e+02   5.7e-02  114.63
    4   5.5e-02  1.1e-03  4.8e-05  1.13e+00   1.893885484e+02   1.893868913e+02   1.1e-03  115.97
    5   1.6e-02  3.1e-04  9.1e-06  1.01e+00   1.894687206e+02   1.894683206e+02   3.1e-04  117.20
    6   1.4e-03  2.7e-05  2.7e-07  1.00e+00   1.895439780e+02   1.895439508e+02   2.7e-05  118.48
    7   1.9e-04  3.8e-06  1.4e-08  1.00e+00   1.895547446e+02   1.895547408e+02   3.8e-06  119.75
    8   1.1e-04  2.1e-06  6.0e-09  1.00e+00   1.895555468e+02   1.895555446e+02   2.1e-06  121.14
    9   7.5e-05  1.5e-06  3.5e-09  1.00e+00   1.895558422e+02   1.895558407e+02   1.5e-06  122.48
    10  2.8e-05  5.6e-07  8.1e-10  1.00e+00   1.895562935e+02   1.895562930e+02   5.6e-07  123.78
    11  4.0e-06  8.0e-08  3.7e-11  1.00e+00   1.895565244e+02   1.895565243e+02   8.0e-08  125.06
    12  2.0e-06  4.0e-08  1.1e-11  1.00e+00   1.895565436e+02   1.895565435e+02   4.0e-08  126.55
    13  1.0e-06  2.0e-08  5.3e-12  1.00e+00   1.895565533e+02   1.895565533e+02   2.0e-08  127.97
    14  5.1e-07  1.0e-08  3.1e-12  1.00e+00   1.895565581e+02   1.895565582e+02   1.0e-08  129.44
    Optimizer terminated. Time: 129.88  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.8955655815e+02    nrm: 4e+02    Viol.  con: 3e-08    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.8955655819e+02    nrm: 1e+02    Viol.  con: 1e-07    var: 2e-10    cones: 0e+00  




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



.. image:: /auto_examples/coil_design/images/sphx_glr_mamba_coil_design_003.png
    :class: sphx-glr-single-img






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 4 minutes  5.728 seconds)


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
