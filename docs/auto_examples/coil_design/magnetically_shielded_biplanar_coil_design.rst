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


    from bfieldtools.mesh_class import MeshWrapper
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
    coil = MeshWrapper(mesh_obj=joined_planes, fix_normals=True)

    # Separate object for shield geometry
    shieldmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/closed_cylinder.stl'), process=True)
    shieldmesh.apply_scale(15)

    shield = MeshWrapper(mesh_obj=shieldmesh, process=True, fix_normals=True)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    face_normals didn't match triangles, ignoring!



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




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_001.png
    :class: sphx-glr-single-img




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

    Computing magnetic field coupling matrix, 3184 vertices by 672 target points... took 0.95 seconds.
    Computing inductance matrix in 2 chunks since 8 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 71.36 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 6930            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 6930            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 6930              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 1.89              dense det. time        : 0.00            
    Factor     - ML order time          : 0.29              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  66.20 
    1   4.0e+01  6.2e-01  2.2e-01  1.05e+00   4.677130957e+01   4.602040536e+01   6.2e-01  66.78 
    2   8.8e+00  1.4e-01  2.0e-02  1.19e+00   7.851531865e+01   7.839683819e+01   1.4e-01  67.34 
    3   4.2e+00  6.6e-02  9.0e-03  1.30e+00   8.060276587e+01   8.055197009e+01   6.6e-02  67.86 
    4   3.2e+00  5.0e-02  7.0e-03  1.01e+00   8.118900965e+01   8.114766773e+01   5.0e-02  68.41 
    5   9.6e-02  1.5e-03  3.8e-05  1.09e+00   8.540388976e+01   8.540282440e+01   1.5e-03  69.12 
    6   1.0e-02  1.6e-04  1.7e-06  1.01e+00   8.545287301e+01   8.545276642e+01   1.6e-04  69.82 
    7   1.1e-03  1.7e-05  5.6e-08  1.00e+00   8.547050333e+01   8.547049238e+01   1.7e-05  70.57 
    8   2.2e-04  3.5e-06  5.3e-09  1.00e+00   8.547222399e+01   8.547222171e+01   3.5e-06  71.33 
    9   1.1e-04  1.7e-06  1.9e-09  1.00e+00   8.547245155e+01   8.547245041e+01   1.7e-06  71.88 
    10  4.2e-06  6.6e-08  1.4e-11  1.00e+00   8.547266623e+01   8.547266619e+01   6.6e-08  72.65 
    11  4.7e-07  7.4e-09  6.7e-13  1.00e+00   8.547267390e+01   8.547267389e+01   7.3e-09  73.60 
    Optimizer terminated. Time: 74.08   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 8.5472673901e+01    nrm: 2e+02    Viol.  con: 3e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 8.5472673895e+01    nrm: 1e+03    Viol.  con: 3e-07    var: 4e-10    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_002.png
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

    Computing scalar potential coupling matrix, 3184 vertices by 962 target points... took 13.00 seconds.
    Computing scalar potential coupling matrix, 962 vertices by 962 target points... took 4.23 seconds.



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





.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_003.png
    :class: sphx-glr-single-img

.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_004.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 962 vertices by 672 target points... took 0.27 seconds.
    This object has no scalar data



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
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 6930            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 6930              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 1.83              dense det. time        : 0.00            
    Factor     - ML order time          : 0.29              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  65.83 
    1   4.0e+01  6.2e-01  3.8e-01  1.08e+00   4.298424044e+01   4.225697886e+01   6.2e-01  66.41 
    2   8.4e+00  1.3e-01  2.0e-02  1.24e+00   6.738470925e+01   6.726360387e+01   1.3e-01  66.97 
    3   1.4e+00  2.2e-02  2.9e-03  1.32e+00   6.906336903e+01   6.904939755e+01   2.2e-02  67.73 
    4   5.8e-01  9.1e-03  7.8e-04  1.06e+00   6.974530407e+01   6.973947021e+01   9.1e-03  68.29 
    5   5.0e-02  7.8e-04  1.8e-05  1.03e+00   7.028043714e+01   7.027992526e+01   7.8e-04  68.97 
    6   2.2e-02  3.4e-04  5.2e-06  9.94e-01   7.031558326e+01   7.031536097e+01   3.4e-04  69.52 
    7   3.2e-04  5.0e-06  9.3e-09  1.00e+00   7.034176335e+01   7.034176003e+01   5.0e-06  70.27 
    8   1.2e-05  1.9e-07  6.9e-11  1.00e+00   7.034214965e+01   7.034214952e+01   1.9e-07  71.02 
    9   4.1e-07  6.3e-09  3.4e-12  1.00e+00   7.034216406e+01   7.034216411e+01   6.3e-09  71.78 
    Optimizer terminated. Time: 72.26   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 7.0342164063e+01    nrm: 1e+02    Viol.  con: 3e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 7.0342164107e+01    nrm: 1e+03    Viol.  con: 5e-07    var: 2e-10    cones: 0e+00  



Plot the newly designed coil windings and field at the target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.j2, N_contours=10)
    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target2 = total_B_coupling @ coil.j2
    mlab.quiver3d(*target_points.T, *B_target2.T)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_005.png
    :class: sphx-glr-single-img




Plot the difference in stream functions


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_data_on_vertices(coil.mesh, np.nan_to_num(100 * (coil.j-coil.j2)/coil.j), figure=f, colorbar=True)

    mlab.colorbar(title='Relative error (%)')



.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_006.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /l/bfieldtools/examples/coil_design/magnetically_shielded_biplanar_coil_design.py:239: RuntimeWarning: invalid value encountered in true_divide
      plot_data_on_vertices(coil.mesh, np.nan_to_num(100 * (coil.j-coil.j2)/coil.j), figure=f, colorbar=True)




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 5 minutes  34.876 seconds)

**Estimated memory usage:**  7854 MB


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
