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





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    SVG path loading unavailable!
    Traceback (most recent call last):
      File "/u/76/zetterr1/unix/.local/lib/python3.6/site-packages/trimesh/path/exchange/svg_io.py", line 18, in <module>
        from svg.path import parse_path
    ModuleNotFoundError: No module named 'svg'



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

    Computing C matrix, 3184 vertices by 512 target points... took 0.84 seconds.



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

    Computing inductance matrix in 2 chunks since 7 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 70.82 seconds.
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
    Factor     - setup time             : 1.78              dense det. time        : 0.00            
    Factor     - ML order time          : 0.37              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.53e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   2.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  73.46 
    1   9.2e+00  3.8e-01  3.3e-01  3.18e-01   8.720442621e+01   8.663738472e+01   3.8e-01  74.14 
    2   6.5e+00  2.7e-01  2.4e-01  3.48e-01   1.226961211e+02   1.222523302e+02   2.7e-01  74.71 
    3   3.4e+00  1.4e-01  1.3e-01  3.68e-01   1.987899646e+02   1.985326339e+02   1.4e-01  75.32 
    4   1.1e+00  4.6e-02  2.1e-02  1.11e+00   3.144686625e+02   3.143947436e+02   4.6e-02  75.90 
    5   4.4e-01  1.8e-02  6.2e-03  9.35e-01   3.445450804e+02   3.445201507e+02   1.8e-02  76.46 
    6   3.2e-02  1.3e-03  1.2e-04  9.84e-01   3.765035939e+02   3.765016483e+02   1.3e-03  77.26 
    7   8.4e-04  3.5e-05  5.6e-07  1.00e+00   3.791913665e+02   3.791913288e+02   3.5e-05  78.02 
    8   3.7e-04  1.5e-05  1.7e-07  1.00e+00   3.792467051e+02   3.792466888e+02   1.5e-05  78.59 
    9   6.2e-06  2.6e-07  4.0e-10  1.00e+00   3.792916856e+02   3.792916853e+02   2.6e-07  79.34 
    10  7.4e-07  3.1e-08  1.7e-11  1.00e+00   3.792923899e+02   3.792923899e+02   3.1e-08  79.92 
    11  3.7e-07  1.5e-08  4.1e-12  1.00e+00   3.792924377e+02   3.792924377e+02   1.5e-08  81.02 
    12  3.9e-07  7.7e-09  2.0e-12  1.00e+00   3.792924617e+02   3.792924617e+02   7.7e-09  81.92 
    13  1.1e-06  3.9e-09  3.8e-12  1.00e+00   3.792924737e+02   3.792924738e+02   3.9e-09  82.91 
    Optimizer terminated. Time: 83.34   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 3.7929247369e+02    nrm: 8e+02    Viol.  con: 2e-08    var: 0e+00    cones: 0e+00  
      Dual.    obj: 3.7929247376e+02    nrm: 4e+02    Viol.  con: 0e+00    var: 2e-10    cones: 0e+00  



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

   **Total running time of the script:** ( 3 minutes  39.183 seconds)

**Estimated memory usage:**  7903 MB


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
