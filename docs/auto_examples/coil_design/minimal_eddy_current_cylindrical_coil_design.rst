.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_minimal_eddy_current_cylindrical_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_minimal_eddy_current_cylindrical_coil_design.py:


Coil with minimal eddy currents
===============================
Compact example of design of a cylindrical coil surrounded by a RF shield, i.e. a conductive surface.
The effects of eddy currents due to inductive interaction with the shield is minimized


.. code-block:: default



    import numpy as np
    from mayavi import mlab
    import trimesh


    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.mesh_properties import mutual_inductance_matrix
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load example coil mesh that is centered on the origin
    coilmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/cylinder.stl'), process=True)

    coilmesh.apply_scale(0.75)

    coilmesh.vertices, coilmesh.faces = trimesh.remesh.subdivide(coilmesh.vertices, coilmesh.faces)

    #Specify offset from origin
    center_offset = np.array([0, 0, 0.75])

    #Apply offset
    coilmesh = trimesh.Trimesh(coilmesh.vertices + center_offset,
                                coilmesh.faces, process=False)

    #Create mesh class object
    coil = MeshWrapper(verts=coilmesh.vertices, tris=coilmesh.faces, fix_normals=True)

    # Separate object for shield geometry
    shield = MeshWrapper(mesh_file=pkg_resources.resource_filename('bfieldtools', 'example_meshes/cylinder.stl'), process=True, fix_normals=True)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    face_normals didn't match triangles, ignoring!
    face_normals didn't match triangles, ignoring!



Set up target  points and plot geometry


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 3])

    sidelength = 0.75 * scaling_factor
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


    #Plot coil, shield and target points

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                    size=(800, 800))

    coil.plot_mesh()
    shield.plot_mesh()
    mlab.points3d(*target_points.T)







.. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_001.png
    :class: sphx-glr-single-img




Compute C matrices that are used to compute the generated magnetic field


.. code-block:: default



    mutual_inductance = mutual_inductance_matrix(coil.mesh, shield.mesh)

    # Take into account the field produced by currents induced into the shield
    # NB! This expression is for instantaneous step-function switching of coil current, see Eq. 18 in G.N. Peeren, 2003.

    shield.coupling = np.linalg.solve(-shield.inductance, mutual_inductance.T)
    secondary_C = shield.B_coupling(target_points) @ shield.coupling





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Estimating 127861 MiB required for 3536 times 904 vertices...
    Computing inductance matrix in 14 chunks since 9563 MiB memory is available...
    Computing potential matrix
    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 32688 MiB required for 904 times 904 vertices...
    Computing inductance matrix in 4 chunks since 9267 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 6.65 seconds.
    Computing magnetic field coupling matrix, 904 vertices by 672 target points... took 0.23 seconds.



Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function
    target_field = np.zeros(target_points.shape)
    target_field[:, 1] = target_field[:, 1] + 1

    target_rel_error = np.zeros_like(target_field)
    target_rel_error[:, 1] += 0.01

    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 1] += 0.001
    target_abs_error[:, 0::2] += 0.005

    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}


    induction_spec = {'coupling':secondary_C, 'abs_error':0.1, 'rel_error':0, 'target':np.zeros(target_field.shape)}





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3536 vertices by 672 target points... took 0.88 seconds.



Run QP solver


.. code-block:: default


    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       [target_spec, induction_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )

    shield.induced_j = shield.coupling @ coil.j






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 500131 MiB required for 3536 times 3536 vertices...
    Computing inductance matrix in 54 chunks since 9289 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 121.30 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 11442           
      Cones                  : 1               
      Scalar variables       : 6755            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 11442           
      Cones                  : 1               
      Scalar variables       : 6755            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 3377
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 11442             conic                  : 3378            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 3.36              dense det. time        : 0.00            
    Factor     - ML order time          : 0.23              GP order time          : 0.00            
    Factor     - nonzeros before factor : 5.70e+06          after factor           : 5.70e+06        
    Factor     - dense dim.             : 0                 flops                  : 9.73e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  472.61
    1   6.1e+01  4.7e-01  1.3e+00  -9.43e-01  1.487988014e+01   1.491559251e+01   4.7e-01  473.69
    2   5.3e+01  4.1e-01  1.2e+00  -8.73e-01  3.514956541e+01   3.542670158e+01   4.1e-01  474.61
    3   1.4e+01  1.1e-01  5.5e-01  -8.55e-01  4.186216296e+02   4.229523709e+02   1.1e-01  475.54
    4   1.0e+01  7.8e-02  4.1e-01  -5.31e-01  9.380575739e+02   9.433819270e+02   7.8e-02  476.48
    5   4.1e+00  3.2e-02  1.7e-01  -3.94e-01  3.451914830e+03   3.458304873e+03   3.2e-02  477.41
    6   8.4e-01  6.6e-03  2.3e-02  8.53e-02   8.645261266e+03   8.647996049e+03   6.6e-03  478.61
    7   5.4e-01  4.2e-03  1.2e-02  1.01e+00   9.127149828e+03   9.128985748e+03   4.2e-03  479.71
    8   4.0e-01  3.1e-03  7.8e-03  8.28e-01   1.001184851e+04   1.001330088e+04   3.1e-03  480.70
    9   1.0e-01  7.9e-04  1.0e-03  8.65e-01   1.204405210e+04   1.204444444e+04   7.9e-04  482.02
    10  6.3e-02  4.9e-04  5.3e-04  9.39e-01   1.238451027e+04   1.238477328e+04   4.9e-04  483.00
    11  2.8e-02  2.2e-04  1.7e-04  8.98e-01   1.275089804e+04   1.275102891e+04   2.2e-04  483.93
    12  2.5e-02  2.0e-04  1.5e-04  8.85e-01   1.278660986e+04   1.278673590e+04   2.0e-04  484.95
    13  5.9e-03  4.6e-05  1.8e-05  8.55e-01   1.308650146e+04   1.308653542e+04   4.6e-05  486.14
    14  1.6e-04  1.3e-06  8.6e-08  9.65e-01   1.319964221e+04   1.319964324e+04   1.3e-06  487.52
    15  1.1e-05  6.9e-07  1.3e-09  9.99e-01   1.320311783e+04   1.320311790e+04   8.2e-08  488.84
    16  8.3e-06  5.5e-07  5.0e-10  1.00e+00   1.320316759e+04   1.320316764e+04   6.5e-08  490.76
    17  6.2e-06  4.0e-07  1.9e-10  1.00e+00   1.320321711e+04   1.320321715e+04   4.8e-08  492.47
    18  2.4e-06  1.6e-07  1.3e-10  1.00e+00   1.320330238e+04   1.320330239e+04   1.9e-08  494.26
    19  2.0e-06  1.3e-07  1.1e-10  1.00e+00   1.320331154e+04   1.320331156e+04   1.6e-08  496.18
    20  1.3e-06  8.6e-08  1.1e-10  1.00e+00   1.320332768e+04   1.320332769e+04   1.0e-08  498.04
    21  2.0e-06  3.3e-08  8.7e-12  1.00e+00   1.320334668e+04   1.320334668e+04   3.9e-09  500.07
    22  1.6e-06  2.5e-08  3.5e-11  1.00e+00   1.320334947e+04   1.320334947e+04   3.0e-09  501.91
    23  5.0e-06  1.5e-08  5.4e-11  1.00e+00   1.320335302e+04   1.320335303e+04   1.9e-09  503.84
    24  4.5e-06  8.4e-09  2.1e-11  1.00e+00   1.320335551e+04   1.320335551e+04   1.0e-09  505.70
    25  8.7e-07  6.1e-09  3.7e-12  1.00e+00   1.320335625e+04   1.320335626e+04   7.7e-10  507.56
    26  1.1e-06  5.4e-09  3.3e-11  1.00e+00   1.320335654e+04   1.320335654e+04   6.8e-10  509.45
    27  1.0e-06  5.3e-09  1.5e-11  1.00e+00   1.320335654e+04   1.320335654e+04   6.8e-10  511.39
    28  1.4e-06  4.9e-09  2.8e-11  1.00e+00   1.320335667e+04   1.320335667e+04   6.3e-10  513.33
    29  1.4e-06  4.5e-09  3.1e-11  1.00e+00   1.320335678e+04   1.320335678e+04   5.9e-10  515.15
    30  1.4e-06  4.5e-09  3.1e-11  1.00e+00   1.320335678e+04   1.320335678e+04   5.9e-10  517.19
    31  1.1e-06  2.3e-09  7.9e-12  1.00e+00   1.320335767e+04   1.320335767e+04   3.0e-10  519.04
    32  1.1e-06  2.3e-09  1.5e-13  1.00e+00   1.320335767e+04   1.320335767e+04   3.0e-10  521.29
    33  5.8e-07  2.1e-09  2.0e-11  1.00e+00   1.320335772e+04   1.320335773e+04   2.8e-10  523.24
    34  1.1e-06  2.1e-09  1.7e-12  1.00e+00   1.320335775e+04   1.320335775e+04   2.7e-10  525.09
    35  1.1e-06  2.1e-09  1.7e-12  1.00e+00   1.320335775e+04   1.320335775e+04   2.7e-10  527.37
    36  1.1e-06  2.1e-09  1.7e-12  1.00e+00   1.320335775e+04   1.320335775e+04   2.7e-10  529.32
    Optimizer terminated. Time: 532.28  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.3203357750e+04    nrm: 3e+04    Viol.  con: 8e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.3203357750e+04    nrm: 5e+05    Viol.  con: 2e-06    var: 2e-08    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default



    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.02)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)


    mlab.title('Coils which minimize the transient effects of conductive shield')





.. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_002.png
    :class: sphx-glr-single-img




For comparison, let's see how the coils look when we ignore the conducting shield


.. code-block:: default



    coil.unshielded_j, coil.unshielded_prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )

    shield.unshielded_induced_j = shield.coupling @ coil.unshielded_j

    loops, loop_values= scalar_contour(coil.mesh, coil.unshielded_j, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.02)

    B_target_unshielded = coil.B_coupling(target_points) @ coil.unshielded_j

    mlab.quiver3d(*target_points.T, *B_target_unshielded.T)

    mlab.title('Coils which ignore the conductive shield')




.. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_003.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 7410            
      Cones                  : 1               
      Scalar variables       : 6755            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 7410            
      Cones                  : 1               
      Scalar variables       : 6755            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 3377
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 7410              conic                  : 3378            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 2.36              dense det. time        : 0.00            
    Factor     - ML order time          : 0.41              GP order time          : 0.00            
    Factor     - nonzeros before factor : 5.70e+06          after factor           : 5.70e+06        
    Factor     - dense dim.             : 0                 flops                  : 7.43e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  138.79
    1   4.0e+01  6.3e-01  2.7e-01  1.58e+00   3.946036018e+01   3.877926361e+01   6.3e-01  139.57
    2   6.5e+00  1.0e-01  1.2e-02  1.13e+00   4.160912452e+01   4.151512039e+01   1.0e-01  140.30
    3   7.9e-01  1.2e-02  9.6e-04  1.41e+00   4.233472661e+01   4.232668490e+01   1.2e-02  141.26
    4   1.3e-01  2.0e-03  6.6e-05  1.06e+00   4.245364929e+01   4.245240209e+01   2.0e-03  142.24
    5   1.9e-02  2.9e-04  3.8e-06  1.01e+00   4.248963943e+01   4.248945954e+01   2.9e-04  143.11
    6   1.6e-03  2.5e-05  9.5e-08  1.00e+00   4.249583959e+01   4.249582414e+01   2.5e-05  144.00
    7   1.2e-03  1.8e-05  5.9e-08  1.00e+00   4.249601118e+01   4.249599990e+01   1.8e-05  144.76
    8   6.5e-04  1.0e-05  2.4e-08  1.00e+00   4.249626489e+01   4.249625868e+01   1.0e-05  145.48
    9   1.7e-05  2.6e-07  1.0e-10  1.00e+00   4.249655242e+01   4.249655225e+01   2.6e-07  146.46
    10  2.9e-07  4.2e-09  2.7e-12  1.00e+00   4.249655983e+01   4.249655988e+01   4.2e-09  147.43
    Optimizer terminated. Time: 148.02  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 4.2496559835e+01    nrm: 9e+01    Viol.  con: 2e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 4.2496559877e+01    nrm: 2e+02    Viol.  con: 3e-07    var: 1e-10    cones: 0e+00  




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 16 minutes  4.760 seconds)

**Estimated memory usage:**  7805 MB


.. _sphx_glr_download_auto_examples_coil_design_minimal_eddy_current_cylindrical_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: minimal_eddy_current_cylindrical_coil_design.py <minimal_eddy_current_cylindrical_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: minimal_eddy_current_cylindrical_coil_design.ipynb <minimal_eddy_current_cylindrical_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
