.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_head_gradient_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_head_gradient_coil_design.py:


Head gradient coil
==================

Example showing a gradient coil designed on the surface of a MEG system helmet


.. code-block:: default



    import numpy as np
    from mayavi import mlab
    import trimesh


    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.magnetic_field_mesh import compute_C
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load simple plane mesh that is centered on the origin
    helmetmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools',
                                                                       'example_meshes/meg_helmet.obj'),
                              process=False)

    #planemesh.apply_scale(scaling_factor)
    #
    ##Specify coil plane geometry
    #center_offset = np.array([0, 0, 0]) * scaling_factor
    #standoff = np.array([0, 4, 0]) * scaling_factor
    #
    ##Create coil plane pairs
    #coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
    #                         planemesh.faces, process=False)
    ##
    #coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
    #                     planemesh.faces, process=False)

    #joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = MeshWrapper(verts=helmetmesh.vertices, tris=helmetmesh.faces, fix_normals=True)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    SVG path loading unavailable!
    Traceback (most recent call last):
      File "/u/76/zetterr1/unix/.local/lib/python3.6/site-packages/trimesh/path/exchange/svg_io.py", line 18, in <module>
        from svg.path import parse_path
    ModuleNotFoundError: No module named 'svg'



Set up target and stray field points.
Here, the target points are on a volumetric grid within a sphere


.. code-block:: default


    center = np.array([0, 0, 0.04]) * scaling_factor

    sidelength = 0.05 * scaling_factor
    n = 12
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








Compute C matrices that are used to compute the generated magnetic field


.. code-block:: default


    coil.C = compute_C(coil.mesh, target_points)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 2044 vertices by 672 target points... took 0.61 seconds.



Specify target field and run solver


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function
    target_field = np.zeros(target_points.shape)
    target_field[:, 0] = target_field[:, 0] + 1 * target_points[:,0]/np.max(target_points[:,0])

    rel_error = np.zeros_like(target_field)
    rel_error[:, 0] += 0.5

    abs_error = np.zeros_like(target_field)
    abs_error[:, 0] += 0.5
    abs_error[:, 1:3] += 0.5


    target_spec = {'C':coil.C, 'rel_error':rel_error, 'abs_error':abs_error, 'target_field':target_field}

    import mosek

    coil.I, prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 1 chunks since 9 GiB memory is available...
    Calculating potentials, chunk 1/1
    Inductance matrix computation took 30.90 seconds.


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 5979            
      Cones                  : 1               
      Scalar variables       : 3893            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 5979            
      Cones                  : 1               
      Scalar variables       : 3893            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 1946
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 5979              conic                  : 1947            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.91              dense det. time        : 0.00            
    Factor     - ML order time          : 0.09              GP order time          : 0.00            
    Factor     - nonzeros before factor : 1.89e+06          after factor           : 1.89e+06        
    Factor     - dense dim.             : 0                 flops                  : 1.75e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.9e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  37.27 
    1   2.2e+01  1.1e-01  6.3e-01  -9.79e-01  6.350121087e+00   1.226904514e+01   1.1e-01  37.61 
    2   6.1e+00  3.2e-02  2.7e-01  -7.93e-01  1.835479781e+01   3.527219599e+01   3.2e-02  37.83 
    3   9.4e-01  4.9e-03  4.8e-02  -4.12e-01  7.638261550e+01   9.977001451e+01   4.9e-03  38.07 
    4   2.8e-02  1.4e-04  3.7e-04  4.57e-01   9.658276130e+01   9.824096869e+01   1.4e-04  38.30 
    5   8.7e-03  4.6e-05  5.9e-05  9.79e-01   7.376772385e+01   7.418242294e+01   4.6e-05  38.54 
    6   6.2e-04  3.3e-06  9.7e-07  9.93e-01   5.956608800e+01   5.958789531e+01   3.3e-06  38.87 
    7   1.6e-04  8.1e-07  1.2e-07  1.00e+00   5.839147102e+01   5.839641393e+01   8.1e-07  39.13 
    8   2.5e-05  1.3e-07  7.1e-09  1.00e+00   5.759453648e+01   5.759522420e+01   1.3e-07  39.39 
    9   7.4e-06  3.8e-08  1.1e-09  1.00e+00   5.749433976e+01   5.749454619e+01   3.8e-08  39.64 
    10  2.5e-06  1.3e-08  2.2e-10  1.00e+00   5.749118663e+01   5.749125846e+01   1.3e-08  39.87 
    11  4.8e-09  2.5e-11  2.8e-14  1.00e+00   5.748651428e+01   5.748651443e+01   2.5e-11  40.20 
    Optimizer terminated. Time: 40.46   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 5.7486514282e+01    nrm: 1e+02    Viol.  con: 3e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 5.7486514427e+01    nrm: 1e+02    Viol.  con: 4e-07    var: 1e-11    cones: 0e+00  



Plot coil windings and magnetic field in target points


.. code-block:: default



    loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=20)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.05/50)

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    mlab.quiver3d(*target_points.T, *B_target.T)

    f.scene.isometric_view()



.. image:: /auto_examples/coil_design/images/sphx_glr_head_gradient_coil_design_001.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 1 minutes  28.528 seconds)

**Estimated memory usage:**  5815 MB


.. _sphx_glr_download_auto_examples_coil_design_head_gradient_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: head_gradient_coil_design.py <head_gradient_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: head_gradient_coil_design.ipynb <head_gradient_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
