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

    Computing magnetic field coupling matrix, 3184 vertices by 672 target points... took 0.78 seconds.
    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 405514 MiB required for 3184 times 3184 vertices...
    Computing inductance matrix in 43 chunks since 9598 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 85.08 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


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
    Factor     - setup time             : 1.66              dense det. time        : 0.00            
    Factor     - ML order time          : 0.16              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  86.72 
    1   4.3e+01  6.7e-01  4.1e-01  7.55e-01   2.942103006e+01   2.861139271e+01   6.7e-01  87.27 
    2   2.6e+01  4.0e-01  1.8e-01  8.10e-01   6.785913679e+01   6.743784167e+01   4.0e-01  87.78 
    3   1.2e+01  1.9e-01  5.8e-02  2.47e+00   8.505386767e+01   8.494298886e+01   1.9e-01  88.29 
    4   8.4e+00  1.3e-01  3.3e-02  1.42e+00   8.907256009e+01   8.900068184e+01   1.3e-01  88.80 
    5   2.6e+00  4.1e-02  5.9e-03  1.28e+00   9.402146526e+01   9.400127813e+01   4.1e-02  89.52 
    6   2.0e+00  3.1e-02  4.1e-03  1.01e+00   9.472179343e+01   9.470586925e+01   3.1e-02  90.03 
    7   1.3e+00  2.0e-02  2.1e-03  8.99e-01   9.590053508e+01   9.589013316e+01   2.0e-02  90.55 
    8   1.9e-01  3.0e-03  1.2e-04  9.71e-01   9.796012024e+01   9.795849811e+01   3.0e-03  91.17 
    9   5.0e-02  7.7e-04  1.5e-05  9.88e-01   9.824073612e+01   9.824029998e+01   7.7e-04  91.79 
    10  4.4e-03  6.9e-05  4.0e-07  9.92e-01   9.835956192e+01   9.835952224e+01   6.9e-05  92.55 
    11  6.0e-04  9.4e-06  2.0e-08  9.98e-01   9.837059727e+01   9.837059185e+01   9.4e-06  93.13 
    12  8.0e-05  1.2e-06  9.6e-10  9.99e-01   9.837225842e+01   9.837225771e+01   1.2e-06  93.73 
    13  2.5e-05  3.9e-07  1.7e-10  1.00e+00   9.837243441e+01   9.837243419e+01   3.9e-07  94.27 
    14  2.6e-06  4.0e-08  4.4e-11  1.00e+00   9.837250705e+01   9.837250684e+01   4.0e-08  94.82 
    15  3.0e-07  4.2e-09  2.0e-13  1.00e+00   9.837251533e+01   9.837251529e+01   2.2e-11  95.53 
    Optimizer terminated. Time: 96.01   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 9.8372515329e+01    nrm: 2e+02    Viol.  con: 8e-12    var: 0e+00    cones: 0e+00  
      Dual.    obj: 9.8372515286e+01    nrm: 1e+03    Viol.  con: 6e-10    var: 2e-10    cones: 0e+00  



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

    Computing scalar potential coupling matrix, 3184 vertices by 962 target points... took 5.57 seconds.
    Computing scalar potential coupling matrix, 962 vertices by 962 target points... took 1.77 seconds.



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

    Computing magnetic field coupling matrix, 962 vertices by 672 target points... took 0.24 seconds.
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
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


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
    Factor     - setup time             : 1.80              dense det. time        : 0.00            
    Factor     - ML order time          : 0.27              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  84.22 
    1   4.3e+01  6.7e-01  2.3e-01  1.08e+00   3.276474332e+01   3.197930167e+01   6.7e-01  84.78 
    2   1.0e+01  1.6e-01  3.2e-02  1.25e+00   7.142638881e+01   7.130956258e+01   1.6e-01  85.35 
    3   4.6e+00  7.2e-02  1.2e-02  1.56e+00   7.589981163e+01   7.585611156e+01   7.2e-02  85.92 
    4   2.7e+00  4.1e-02  5.1e-03  1.22e+00   7.831447027e+01   7.829093435e+01   4.1e-02  86.44 
    5   4.9e-01  7.5e-03  4.3e-04  1.13e+00   8.025516651e+01   8.025123784e+01   7.5e-03  87.16 
    6   2.0e-01  3.1e-03  1.2e-04  1.02e+00   8.063674999e+01   8.063511283e+01   3.1e-03  87.68 
    7   3.6e-02  5.5e-04  8.3e-06  9.93e-01   8.087222785e+01   8.087193264e+01   5.5e-04  88.31 
    8   2.1e-02  3.3e-04  3.8e-06  1.00e+00   8.089438000e+01   8.089420499e+01   3.3e-04  88.83 
    9   1.8e-02  2.7e-04  2.9e-06  9.91e-01   8.089997565e+01   8.089983003e+01   2.7e-04  89.35 
    10  2.3e-03  3.6e-05  1.4e-07  1.00e+00   8.092364663e+01   8.092362743e+01   3.6e-05  89.89 
    11  9.6e-06  1.5e-07  3.2e-11  1.00e+00   8.092736660e+01   8.092736653e+01   1.5e-07  90.61 
    12  1.6e-06  2.5e-08  9.0e-13  1.00e+00   8.092737936e+01   8.092737936e+01   2.5e-08  91.13 
    Optimizer terminated. Time: 91.61   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 8.0927379364e+01    nrm: 2e+02    Viol.  con: 8e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 8.0927379361e+01    nrm: 1e+03    Viol.  con: 7e-07    var: 1e-10    cones: 0e+00  



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

   **Total running time of the script:** ( 6 minutes  18.161 seconds)

**Estimated memory usage:**  4093 MB


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
