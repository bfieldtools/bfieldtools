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
    from bfieldtools.magnetic_field_mesh import compute_C
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix
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


    coil.C = compute_C(coil.mesh, target_points)
    shield.C = compute_C(shield.mesh, target_points)

    mutual_inductance = mutual_inductance_matrix(coil.mesh, shield.mesh)

    # Take into account the field produced by currents induced into the shield
    # NB! This expression is for instantaneous step-function switching of coil current, see Eq. 18 in G.N. Peeren, 2003.

    shield.coupling = np.linalg.solve(-shield.inductance, mutual_inductance.T)
    secondary_C = (shield.C.transpose((0,2,1)) @ shield.coupling).transpose((0,2,1))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 3536 vertices by 672 target points... took 1.03 seconds.
    Computing C matrix, 904 vertices by 672 target points... took 0.23 seconds.
    Calculating potentials
    Inserting stuff into M-matrix
    Computing inductance matrix in 1 chunks since 8 GiB memory is available...
    Calculating potentials, chunk 1/1
    Inductance matrix computation took 5.59 seconds.



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

    target_spec = {'C':coil.C, 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target_field':target_field}


    induction_spec = {'C':secondary_C, 'abs_error':0.1, 'rel_error':0, 'target_field':np.zeros(target_field.shape)}







Run QP solver


.. code-block:: default


    import mosek

    coil.I, prob = optimize_streamfunctions(coil,
                                       [target_spec, induction_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )

    shield.induced_I = shield.coupling @ coil.I






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 2 chunks since 8 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 90.48 seconds.


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
    Factor     - setup time             : 3.54              dense det. time        : 0.00            
    Factor     - ML order time          : 0.42              GP order time          : 0.00            
    Factor     - nonzeros before factor : 5.70e+06          after factor           : 5.70e+06        
    Factor     - dense dim.             : 0                 flops                  : 9.73e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  419.42
    1   5.7e+01  4.5e-01  1.3e+00  -9.57e-01  8.859804290e+00   9.030546636e+00   4.5e-01  420.54
    2   4.0e+01  3.1e-01  1.0e+00  -8.86e-01  3.356983293e+01   3.452221098e+01   3.1e-01  421.63
    3   2.4e+01  1.9e-01  7.6e-01  -8.34e-01  1.181609808e+02   1.205028447e+02   1.9e-01  422.78
    4   2.0e+01  1.6e-01  6.8e-01  -7.36e-01  1.578027546e+02   1.606427513e+02   1.6e-01  423.79
    5   1.1e+01  8.8e-02  4.4e-01  -6.92e-01  6.952681970e+02   7.000013282e+02   8.8e-02  424.80
    6   7.8e+00  6.1e-02  3.2e-01  -4.96e-01  1.027360428e+03   1.033096344e+03   6.1e-02  425.77
    7   3.9e+00  3.1e-02  1.6e-01  -3.22e-01  3.350722063e+03   3.356814020e+03   3.1e-02  426.71
    8   6.5e-01  5.1e-03  1.8e-02  9.43e-02   7.789314649e+03   7.792144350e+03   5.1e-03  428.05
    9   5.4e-01  4.2e-03  1.4e-02  9.98e-01   7.910329862e+03   7.912760121e+03   4.2e-03  429.00
    10  4.6e-01  3.6e-03  1.1e-02  8.88e-01   8.261351733e+03   8.263517198e+03   3.6e-03  429.93
    11  1.1e-01  8.7e-04  1.2e-03  9.19e-01   1.005073104e+04   1.005118914e+04   8.7e-04  431.21
    12  6.5e-02  5.0e-04  5.6e-04  9.66e-01   1.045041260e+04   1.045069816e+04   5.0e-04  432.13
    13  3.6e-02  2.8e-04  2.4e-04  9.27e-01   1.070048512e+04   1.070065546e+04   2.8e-04  433.12
    14  2.1e-02  1.7e-04  1.2e-04  9.42e-01   1.083380734e+04   1.083391894e+04   1.7e-04  434.15
    15  6.2e-03  4.9e-05  1.9e-05  8.49e-01   1.104262837e+04   1.104266376e+04   4.9e-05  435.15
    16  1.5e-04  1.2e-06  7.6e-08  9.62e-01   1.113623823e+04   1.113623916e+04   1.2e-06  436.43
    17  7.5e-06  2.3e-07  7.9e-10  9.99e-01   1.113875172e+04   1.113875177e+04   5.8e-08  437.70
    18  2.0e-07  1.4e-07  6.3e-12  1.00e+00   1.113888085e+04   1.113888085e+04   9.7e-10  438.97
    19  2.0e-07  1.4e-07  6.3e-12  1.00e+00   1.113888085e+04   1.113888085e+04   9.7e-10  441.01
    20  7.7e-07  1.3e-07  8.9e-13  1.00e+00   1.113888081e+04   1.113888081e+04   8.6e-10  442.83
    21  2.3e-06  7.4e-08  3.3e-12  1.00e+00   1.113888078e+04   1.113888078e+04   5.0e-10  444.58
    22  2.3e-06  7.4e-08  3.3e-12  1.00e+00   1.113888078e+04   1.113888078e+04   5.0e-10  446.52
    23  1.8e-06  6.4e-08  2.9e-12  1.00e+00   1.113888118e+04   1.113888118e+04   4.3e-10  448.30
    24  2.6e-06  1.6e-08  2.0e-12  1.00e+00   1.113888162e+04   1.113888162e+04   1.0e-10  450.10
    25  2.6e-06  1.6e-08  2.0e-12  1.00e+00   1.113888162e+04   1.113888162e+04   1.0e-10  452.36
    26  2.6e-06  1.5e-08  4.4e-12  1.00e+00   1.113888164e+04   1.113888164e+04   9.4e-11  454.26
    27  2.6e-06  1.5e-08  4.4e-12  1.00e+00   1.113888164e+04   1.113888164e+04   9.4e-11  456.47
    28  3.0e-06  1.4e-08  4.7e-12  1.00e+00   1.113888166e+04   1.113888166e+04   8.3e-11  458.32
    29  3.0e-06  1.3e-08  6.0e-12  1.00e+00   1.113888167e+04   1.113888168e+04   7.8e-11  460.23
    30  3.0e-06  1.3e-08  6.0e-12  1.00e+00   1.113888167e+04   1.113888168e+04   7.8e-11  462.50
    Optimizer terminated. Time: 465.40  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.1138881674e+04    nrm: 2e+04    Viol.  con: 2e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.1138881675e+04    nrm: 4e+05    Viol.  con: 2e-06    var: 9e-08    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default



    loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.02)

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    mlab.quiver3d(*target_points.T, *B_target.T)


    mlab.title('Coils which minimize the transient effects of conductive shield')





.. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_002.png
    :class: sphx-glr-single-img




For comparison, let's see how the coils look when we ignore the conducting shield


.. code-block:: default



    coil.unshielded_I, coil.unshielded_prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )

    shield.unshielded_induced_I = shield.coupling @ coil.unshielded_I

    loops, loop_values= scalar_contour(coil.mesh, coil.unshielded_I, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.02)

    B_target_unshielded = coil.C.transpose([0, 2, 1]) @ coil.unshielded_I

    mlab.quiver3d(*target_points.T, *B_target_unshielded.T)

    mlab.title('Coils which ignore the conductive shield')




.. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_003.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none



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
    Factor     - setup time             : 2.37              dense det. time        : 0.00            
    Factor     - ML order time          : 0.41              GP order time          : 0.00            
    Factor     - nonzeros before factor : 5.70e+06          after factor           : 5.70e+06        
    Factor     - dense dim.             : 0                 flops                  : 7.43e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  106.00
    1   3.8e+01  6.0e-01  2.3e-01  1.63e+00   3.716704927e+01   3.652271215e+01   6.0e-01  106.77
    2   5.9e+00  9.1e-02  1.1e-02  1.14e+00   3.467322805e+01   3.458895984e+01   9.1e-02  107.48
    3   1.0e+00  1.6e-02  1.1e-03  1.40e+00   3.557447001e+01   3.556329306e+01   1.6e-02  108.29
    4   1.4e-01  2.2e-03  8.1e-05  1.08e+00   3.569632267e+01   3.569495849e+01   2.2e-03  109.13
    5   6.7e-02  1.0e-03  2.6e-05  1.01e+00   3.571818972e+01   3.571755284e+01   1.0e-03  109.87
    6   1.1e-02  1.6e-04  1.6e-06  1.01e+00   3.571619267e+01   3.571609181e+01   1.6e-04  110.75
    7   1.6e-03  2.6e-05  9.8e-08  1.00e+00   3.571984480e+01   3.571982909e+01   2.6e-05  111.55
    8   7.3e-05  1.1e-06  2.1e-09  1.00e+00   3.572043680e+01   3.572043496e+01   1.1e-06  112.51
    9   5.0e-05  7.8e-07  1.2e-09  1.00e+00   3.572044678e+01   3.572044552e+01   7.8e-07  113.25
    10  2.6e-05  4.0e-07  4.5e-10  1.00e+00   3.572045651e+01   3.572045586e+01   4.0e-07  113.98
    11  4.7e-07  7.4e-09  4.9e-12  1.00e+00   3.572046686e+01   3.572046692e+01   7.4e-09  114.88
    Optimizer terminated. Time: 115.46  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 3.5720466860e+01    nrm: 7e+01    Viol.  con: 3e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 3.5720466918e+01    nrm: 1e+02    Viol.  con: 5e-07    var: 2e-10    cones: 0e+00  




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 13 minutes  30.264 seconds)

**Estimated memory usage:**  10256 MB


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
