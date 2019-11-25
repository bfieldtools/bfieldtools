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
    from bfieldtools.mesh_inductance import mutual_inductance_matrix
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

    Calculating potentials
    Inserting stuff into M-matrix
    Computing inductance matrix in 1 chunks since 7 GiB memory is available...
    Calculating potentials, chunk 1/1
    Inductance matrix computation took 6.00 seconds.
    Computing C matrix, 904 vertices by 672 target points... took 0.26 seconds.



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

    Computing C matrix, 3536 vertices by 672 target points... took 1.16 seconds.



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

    Computing inductance matrix in 3 chunks since 6 GiB memory is available...
    Calculating potentials, chunk 1/3
    Calculating potentials, chunk 2/3
    Calculating potentials, chunk 3/3
    Inductance matrix computation took 92.89 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


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
    Factor     - setup time             : 3.56              dense det. time        : 0.00            
    Factor     - ML order time          : 0.42              GP order time          : 0.00            
    Factor     - nonzeros before factor : 5.70e+06          after factor           : 5.70e+06        
    Factor     - dense dim.             : 0                 flops                  : 9.73e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  420.80
    1   5.7e+01  4.5e-01  1.3e+00  -9.57e-01  8.859804290e+00   9.030546622e+00   4.5e-01  421.91
    2   4.0e+01  3.1e-01  1.0e+00  -8.86e-01  3.356983293e+01   3.452221096e+01   3.1e-01  422.86
    3   2.4e+01  1.9e-01  7.6e-01  -8.34e-01  1.181609808e+02   1.205028447e+02   1.9e-01  423.78
    4   2.0e+01  1.6e-01  6.8e-01  -7.36e-01  1.578027546e+02   1.606427513e+02   1.6e-01  424.72
    5   1.1e+01  8.8e-02  4.4e-01  -6.92e-01  6.952681970e+02   7.000013282e+02   8.8e-02  425.66
    6   7.8e+00  6.1e-02  3.2e-01  -4.96e-01  1.027360428e+03   1.033096344e+03   6.1e-02  426.62
    7   3.9e+00  3.1e-02  1.6e-01  -3.22e-01  3.350722062e+03   3.356814020e+03   3.1e-02  427.56
    8   6.5e-01  5.1e-03  1.8e-02  9.43e-02   7.789314650e+03   7.792144350e+03   5.1e-03  428.87
    9   5.4e-01  4.2e-03  1.4e-02  9.98e-01   7.910329862e+03   7.912760121e+03   4.2e-03  429.82
    10  4.6e-01  3.6e-03  1.1e-02  8.88e-01   8.261351733e+03   8.263517198e+03   3.6e-03  430.75
    11  1.1e-01  8.7e-04  1.2e-03  9.19e-01   1.005073104e+04   1.005118914e+04   8.7e-04  432.05
    12  6.5e-02  5.0e-04  5.6e-04  9.66e-01   1.045041260e+04   1.045069816e+04   5.0e-04  433.03
    13  3.6e-02  2.8e-04  2.4e-04  9.27e-01   1.070048511e+04   1.070065545e+04   2.8e-04  433.97
    14  2.1e-02  1.7e-04  1.2e-04  9.42e-01   1.083380734e+04   1.083391894e+04   1.7e-04  434.92
    15  6.2e-03  4.9e-05  1.9e-05  8.49e-01   1.104262837e+04   1.104266375e+04   4.9e-05  435.93
    16  1.5e-04  1.2e-06  7.6e-08  9.62e-01   1.113623823e+04   1.113623916e+04   1.2e-06  437.26
    17  7.5e-06  2.4e-07  8.2e-10  9.99e-01   1.113875172e+04   1.113875177e+04   5.8e-08  438.56
    18  1.2e-07  5.4e-08  1.6e-11  1.00e+00   1.113888085e+04   1.113888085e+04   9.7e-10  439.87
    19  2.0e-07  4.1e-08  5.3e-12  1.00e+00   1.113888077e+04   1.113888077e+04   7.4e-10  441.76
    20  3.5e-07  2.7e-08  8.1e-13  1.00e+00   1.113888081e+04   1.113888081e+04   4.9e-10  443.53
    21  1.6e-07  2.5e-08  2.6e-12  1.00e+00   1.113888106e+04   1.113888106e+04   4.5e-10  445.31
    22  1.6e-06  1.5e-08  1.1e-12  1.00e+00   1.113888137e+04   1.113888137e+04   2.6e-10  447.11
    23  7.8e-07  1.1e-08  6.0e-12  1.00e+00   1.113888159e+04   1.113888159e+04   1.4e-10  448.96
    24  7.8e-07  1.1e-08  6.0e-12  1.00e+00   1.113888159e+04   1.113888159e+04   1.4e-10  450.89
    25  8.2e-07  1.1e-08  6.9e-14  1.00e+00   1.113888159e+04   1.113888159e+04   1.4e-10  452.87
    26  8.2e-07  1.1e-08  6.9e-14  1.00e+00   1.113888159e+04   1.113888159e+04   1.4e-10  454.91
    27  7.7e-07  1.1e-08  7.7e-12  1.00e+00   1.113888159e+04   1.113888159e+04   1.4e-10  456.85
    28  7.7e-07  1.1e-08  7.7e-12  1.00e+00   1.113888159e+04   1.113888159e+04   1.4e-10  459.02
    29  1.1e-07  8.7e-09  4.9e-12  1.00e+00   1.113888165e+04   1.113888165e+04   1.0e-10  461.04
    30  1.1e-07  8.7e-09  4.9e-12  1.00e+00   1.113888165e+04   1.113888165e+04   1.0e-10  463.08
    31  1.8e-06  7.9e-09  1.4e-12  1.00e+00   1.113888170e+04   1.113888170e+04   7.9e-11  465.02
    32  1.8e-06  7.9e-09  1.4e-12  1.00e+00   1.113888170e+04   1.113888170e+04   7.9e-11  467.32
    33  1.8e-06  7.8e-09  6.0e-12  1.00e+00   1.113888170e+04   1.113888170e+04   7.9e-11  469.47
    34  1.8e-06  7.6e-09  8.7e-12  1.00e+00   1.113888170e+04   1.113888170e+04   7.9e-11  471.51
    35  1.8e-06  7.6e-09  2.5e-12  1.00e+00   1.113888170e+04   1.113888170e+04   7.9e-11  473.54
    36  1.9e-06  7.5e-09  9.3e-13  1.00e+00   1.113888170e+04   1.113888170e+04   7.9e-11  475.64
    Optimizer terminated. Time: 479.01  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.1138881700e+04    nrm: 2e+04    Viol.  con: 2e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.1138881700e+04    nrm: 4e+05    Viol.  con: 2e-06    var: 1e-07    cones: 0e+00  



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
    Factor     - setup time             : 2.52              dense det. time        : 0.00            
    Factor     - ML order time          : 0.47              GP order time          : 0.00            
    Factor     - nonzeros before factor : 5.70e+06          after factor           : 5.70e+06        
    Factor     - dense dim.             : 0                 flops                  : 7.43e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  104.86
    1   3.8e+01  6.0e-01  2.3e-01  1.63e+00   3.716704927e+01   3.652271215e+01   6.0e-01  105.65
    2   5.9e+00  9.1e-02  1.1e-02  1.14e+00   3.467322805e+01   3.458895984e+01   9.1e-02  106.39
    3   1.0e+00  1.6e-02  1.1e-03  1.40e+00   3.557447001e+01   3.556329306e+01   1.6e-02  107.21
    4   1.4e-01  2.2e-03  8.1e-05  1.08e+00   3.569632267e+01   3.569495849e+01   2.2e-03  108.11
    5   6.7e-02  1.0e-03  2.6e-05  1.01e+00   3.571818972e+01   3.571755284e+01   1.0e-03  108.98
    6   1.1e-02  1.6e-04  1.6e-06  1.01e+00   3.571619267e+01   3.571609181e+01   1.6e-04  110.10
    7   1.6e-03  2.6e-05  9.8e-08  1.00e+00   3.571984480e+01   3.571982910e+01   2.6e-05  111.14
    8   7.3e-05  1.1e-06  8.4e-10  1.00e+00   3.572043680e+01   3.572043618e+01   1.1e-06  112.17
    9   5.0e-05  7.8e-07  4.7e-10  1.00e+00   3.572044677e+01   3.572044635e+01   7.8e-07  112.93
    10  2.6e-05  4.0e-07  1.8e-10  1.00e+00   3.572045650e+01   3.572045628e+01   4.0e-07  113.70
    11  4.7e-07  7.4e-09  1.7e-12  1.00e+00   3.572046686e+01   3.572046684e+01   7.4e-09  114.63
    Optimizer terminated. Time: 115.23  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 3.5720466860e+01    nrm: 7e+01    Viol.  con: 3e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 3.5720466841e+01    nrm: 1e+02    Viol.  con: 5e-07    var: 7e-11    cones: 0e+00  




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 14 minutes  26.880 seconds)

**Estimated memory usage:**  7294 MB


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
