.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_magnetically_shielded_biplanar_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_magnetically_shielded_biplanar_coil_design.py:


Magnetically shielded  coil
===========================
Compact example of design of a biplanar coil within a cylindrical shield.
The effect of the shield is prospectively taken into account while designing the coil.
The coil is positioned close to the end of the shield to demonstrate the effect


.. code-block:: default



    import numpy as np
    from mayavi import mlab
    import trimesh


    from bfieldtools.mesh_class import Conductor
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops, plot_data_on_vertices

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

    planemesh.apply_scale(scaling_factor)

    #Specify coil plane geometry
    center_offset = np.array([9, 0, 0]) * scaling_factor
    standoff = np.array([0, 4, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                             planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                         planemesh.faces, process=False)

    joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = Conductor(mesh_obj=joined_planes, fix_normals=True)

    # Separate object for shield geometry
    shieldmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/closed_cylinder.stl'), process=True)
    shieldmesh.apply_scale(15)

    shield = Conductor(mesh_obj=shieldmesh, process=True, fix_normals=True)








Set up target  points and plot geometry


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere
    # Set up target and stray field points

    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([9, 0, 0]) * scaling_factor

    sidelength = 3 * scaling_factor
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

    coil.plot_mesh(representation='surface')
    shield.plot_mesh()
    mlab.points3d(*target_points.T)




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_001.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_003.png
            :class: sphx-glr-multi-img




Let's design a coil without taking the magnetic shield into account


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function
    target_field = np.zeros(target_points.shape)
    target_field[:, 1] = target_field[:, 1] + 1 # Homogeneous Z-field


    target_rel_error = np.zeros_like(target_field)
    target_rel_error[:, 0] += 0.01

    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 0] += 0.001
    target_abs_error[:, 1:3] += 0.005

    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}

    import mosek

    coil.j, coil.prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 672 target points... took 0.93 seconds.
    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 100 chunks (7772 MiB memory free), when approx_far=True using more chunks is faster...
    Computing potential matrix
    Inductance matrix computation took 50.77 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 6930            
      Cones                  : 1               
      Scalar variables       : 6083            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 6930            
      Cones                  : 1               
      Scalar variables       : 6083            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 6930              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 1.82              dense det. time        : 0.00            
    Factor     - ML order time          : 0.29              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  78.38 
    1   3.9e+01  6.0e-01  2.2e-01  9.32e-01   5.048078588e+01   4.972913092e+01   6.0e-01  78.92 
    2   9.2e+00  1.4e-01  2.7e-02  1.15e+00   8.981665197e+01   8.968983538e+01   1.4e-01  79.45 
    3   4.3e+00  6.7e-02  1.1e-02  1.30e+00   9.209533675e+01   9.204271073e+01   6.7e-02  79.98 
    4   1.3e+00  2.0e-02  2.3e-03  1.05e+00   9.637250546e+01   9.635715466e+01   2.0e-02  80.60 
    5   1.3e-01  2.0e-03  7.5e-05  1.07e+00   9.809905145e+01   9.809772029e+01   2.0e-03  81.27 
    6   1.3e-02  2.1e-04  2.6e-06  1.00e+00   9.829109041e+01   9.829095050e+01   2.1e-04  81.90 
    7   7.2e-03  1.1e-04  1.0e-06  9.99e-01   9.830510384e+01   9.830502787e+01   1.1e-04  82.42 
    8   9.8e-04  1.5e-05  5.2e-08  1.00e+00   9.831974067e+01   9.831973029e+01   1.5e-05  83.07 
    9   4.3e-04  6.7e-06  1.5e-08  1.00e+00   9.832108213e+01   9.832107763e+01   6.7e-06  83.57 
    10  2.6e-05  4.0e-07  1.5e-10  1.00e+00   9.832206844e+01   9.832206828e+01   4.0e-07  84.16 
    11  9.1e-06  4.6e-08  1.9e-11  1.00e+00   9.832212605e+01   9.832212597e+01   4.6e-08  84.84 
    12  5.8e-06  3.0e-08  2.0e-11  1.00e+00   9.832212855e+01   9.832212843e+01   3.0e-08  85.36 
    13  2.8e-06  1.4e-08  3.0e-11  1.00e+00   9.832213110e+01   9.832213085e+01   1.4e-08  85.85 
    Optimizer terminated. Time: 86.33   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 9.8322131100e+01    nrm: 2e+02    Viol.  con: 6e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 9.8322130854e+01    nrm: 1e+03    Viol.  con: 5e-07    var: 2e-09    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_004.png
    :class: sphx-glr-single-img




Now, let's compute the effect of the shield on the field produced by the coil


.. code-block:: default


    # Points slightly inside the shield
    d = np.mean(np.diff(shield.mesh.vertices[shield.mesh.faces[:,0:2]],axis=1), axis=0)/10
    points = shield.mesh.vertices - d*shield.mesh.vertex_normals

    # Calculate primary potential matrix at the shield surface
    P_prim = coil.U_coupling(points)

    # Calculate linear collocation BEM matrix
    P_bem = shield.U_coupling(points)

    # Recalculate diag elements according to de Munck paper
    #for diag_index in range(P_bem.shape[0]):
    #    P_bem[diag_index, diag_index] = 0
    #    P_bem[diag_index, diag_index] = -P_bem[diag_index, :].sum()

    # Matrix misses one rank, make it invertible
    # by rank-one update (sets potential of constant dipole layer)
    #P_bem += np.ones(P_bem.shape)/P_bem.shape[0]


    # Solve equivalent stream function for the perfect linear mu-metal layer.
    # This is the equivalent surface current in the shield that would cause its
    # scalar magnetic potential to be constant
    shield.j =  np.linalg.solve(P_bem, P_prim @ coil.j)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing scalar potential coupling matrix, 3184 vertices by 962 target points... took 5.73 seconds.
    Computing scalar potential coupling matrix, 962 vertices by 962 target points... took 1.78 seconds.



Plot the difference in field when taking the shield into account


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    B_target = coil.B_coupling(target_points) @ coil.j

    B_target_w_shield = coil.B_coupling(target_points) @ coil.j + shield.B_coupling(target_points) @ shield.j

    B_quiver = mlab.quiver3d(*target_points.T, *(B_target_w_shield - B_target).T, colormap='viridis', mode='arrow')
    f.scene.isometric_view()
    mlab.colorbar(B_quiver, title='Difference in magnetic field (a.u.)')

    import seaborn as sns
    import matplotlib.pyplot as plt




    fig, axes = plt.subplots(1, 3, figsize=(10, 4))

    fig.suptitle('Component-wise effect of magnetic shield on target field amplitude distribution')
    for ax_idx, ax in enumerate(axes):

        sns.distplot(B_target[:, ax_idx], label='Without shield', ax=ax)
        sns.distplot(B_target_w_shield[:, ax_idx], label='With shield', ax=ax)
        ax.set_xlabel('Magnetic field (a.u.)')

        if ax_idx == 2:
            ax.legend()

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])





.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_005.png
    :class: sphx-glr-single-img

.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_006.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 962 vertices by 672 target points... took 0.27 seconds.
    This object has no scalar data
    /u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/ipykernel/kernelbase.py:19: VisibleDeprecationWarning: zmq.eventloop.minitornado is deprecated in pyzmq 14.0 and will be removed.
        Install tornado itself to use zmq with the tornado IOLoop.
    
      from jupyter_client.session import utcnow as now



Let's redesign the coil taking the shield into account prospectively


.. code-block:: default


    shield.coupling = np.linalg.solve(P_bem, P_prim)

    secondary_B_coupling = shield.B_coupling(target_points) @ shield.coupling

    total_B_coupling = coil.B_coupling(target_points) + secondary_B_coupling

    target_spec_w_shield = {'coupling':total_B_coupling, 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}


    coil.j2, coil.prob2 = optimize_streamfunctions(coil,
                                       [target_spec_w_shield],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )





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
      Constraints            : 6930            
      Cones                  : 1               
      Scalar variables       : 6083            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 6930            
      Cones                  : 1               
      Scalar variables       : 6083            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 6930              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 1.76              dense det. time        : 0.00            
    Factor     - ML order time          : 0.27              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  90.09 
    1   3.7e+01  5.7e-01  3.8e-01  8.88e-01   4.754179191e+01   4.682655503e+01   5.7e-01  90.62 
    2   8.1e+00  1.3e-01  2.1e-02  1.20e+00   7.723939039e+01   7.711370701e+01   1.3e-01  91.16 
    3   1.3e+00  2.0e-02  2.9e-03  1.32e+00   7.934876951e+01   7.933475709e+01   2.0e-02  91.88 
    4   6.0e-01  9.3e-03  8.7e-04  1.06e+00   8.012501237e+01   8.011870723e+01   9.3e-03  92.40 
    5   4.3e-02  6.7e-04  1.6e-05  1.03e+00   8.082162749e+01   8.082116350e+01   6.7e-04  93.09 
    6   3.5e-02  5.5e-04  1.2e-05  9.84e-01   8.083297867e+01   8.083259415e+01   5.5e-04  93.61 
    7   1.1e-02  1.7e-04  1.9e-06  1.00e+00   8.086988981e+01   8.086977400e+01   1.7e-04  94.14 
    8   1.3e-04  2.0e-06  2.5e-09  1.00e+00   8.088581305e+01   8.088581166e+01   2.0e-06  94.72 
    9   4.1e-06  6.4e-08  1.2e-11  1.00e+00   8.088600442e+01   8.088600439e+01   6.3e-08  95.43 
    10  1.0e-06  1.6e-08  3.6e-12  1.00e+00   8.088600933e+01   8.088600930e+01   1.6e-08  95.93 
    Optimizer terminated. Time: 96.39   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 8.0886009329e+01    nrm: 2e+02    Viol.  con: 7e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 8.0886009305e+01    nrm: 1e+03    Viol.  con: 1e-06    var: 1e-10    cones: 0e+00  



Plot the newly designed coil windings and field at the target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.j2, N_contours=10)
    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target2 = total_B_coupling @ coil.j2
    mlab.quiver3d(*target_points.T, *B_target2.T)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_007.png
    :class: sphx-glr-single-img




Plot the difference in stream functions


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_data_on_vertices(coil.mesh, np.nan_to_num(100 * (coil.j-coil.j2)/coil.j), figure=f, colorbar=True)

    mlab.colorbar(title='Relative error (%)')



.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_008.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools/examples/coil_design/magnetically_shielded_biplanar_coil_design.py:239: RuntimeWarning: invalid value encountered in true_divide
      plot_data_on_vertices(coil.mesh, np.nan_to_num(100 * (coil.j-coil.j2)/coil.j), figure=f, colorbar=True)




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 6 minutes  40.532 seconds)


.. _sphx_glr_download_auto_examples_coil_design_magnetically_shielded_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: magnetically_shielded_biplanar_coil_design.py <magnetically_shielded_biplanar_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: magnetically_shielded_biplanar_coil_design.ipynb <magnetically_shielded_biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
