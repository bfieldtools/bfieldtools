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


    from bfieldtools.mesh_class import Conductor, StreamFunction
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
    coil = Conductor(mesh_obj=joined_planes,
                     fix_normals=True,
                     basis_name='inner')

    # Separate object for shield geometry
    shieldmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/closed_cylinder_remeshed.stl'), process=True)
    shieldmesh.apply_scale(15)

    shield = Conductor(mesh_obj=shieldmesh,
                       process=True,
                       fix_normals=True,
                       basis_name='vertex')








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
    shield.plot_mesh(representation='surface', cull_front=True, color=(0.9, 0.9, 0.9))
    mlab.points3d(*target_points.T)


    f.scene.isometric_view()
    f.scene.camera.zoom(1.2)





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
    target_field[:, 0] = target_field[:, 0] + 1 # Homogeneous Y-field



    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 0] += 0.005
    target_abs_error[:, 1:3] += 0.01

    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':0, 'abs_error':target_abs_error, 'target':target_field}

    import mosek

    coil.s, coil.prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 672 target points... took 0.78 seconds.
    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 80 chunks (11621 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 42.43 seconds.
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
    Factor     - setup time             : 1.82              dense det. time        : 0.00            
    Factor     - ML order time          : 0.28              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  85.10 
    1   6.2e+01  4.9e-01  8.3e-01  -1.88e-01  1.238908380e+02   1.233494816e+02   4.9e-01  85.65 
    2   2.8e+01  2.2e-01  2.7e-01  -7.58e-02  4.723125710e+02   4.720214058e+02   2.2e-01  86.20 
    3   6.6e+00  5.1e-02  3.4e-02  1.16e+00   7.663818856e+02   7.663394148e+02   5.1e-02  86.91 
    4   1.2e+00  9.7e-03  2.7e-03  1.01e+00   8.627010809e+02   8.626916977e+02   9.7e-03  87.63 
    5   1.7e-01  1.3e-03  1.4e-04  9.78e-01   8.843341332e+02   8.843330963e+02   1.3e-03  88.18 
    6   2.1e-02  1.6e-04  6.4e-06  1.00e+00   8.876119122e+02   8.876117890e+02   1.6e-04  88.81 
    7   2.6e-03  2.0e-05  2.7e-07  1.00e+00   8.880900344e+02   8.880900201e+02   2.0e-05  89.35 
    8   3.5e-04  2.7e-06  1.3e-08  1.00e+00   8.881465514e+02   8.881465495e+02   2.7e-06  89.86 
    9   2.8e-06  1.7e-08  6.2e-11  1.00e+00   8.881555063e+02   8.881555084e+02   1.1e-09  90.59 
    10  1.4e-05  7.8e-09  1.5e-11  1.00e+00   8.881555079e+02   8.881555087e+02   5.3e-10  91.62 
    Optimizer terminated. Time: 92.09   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 8.8815550795e+02    nrm: 2e+03    Viol.  con: 8e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 8.8815550866e+02    nrm: 6e+03    Viol.  con: 2e-08    var: 3e-09    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.s.vert, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.s

    mlab.quiver3d(*target_points.T, *B_target.T, mode='arrow', scale_factor=0.75)




    f.scene.isometric_view()
    f.scene.camera.zoom(0.95)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_004.png
    :class: sphx-glr-single-img




Now, let's compute the effect of the shield on the field produced by the coil


.. code-block:: default


    # Points slightly inside the shield
    d = np.mean(np.diff(shield.mesh.vertices[shield.mesh.faces[:,0:2]],axis=1), axis=0)/10
    points = shield.mesh.vertices - d*shield.mesh.vertex_normals


    # Solve equivalent stream function for the perfect linear mu-metal layer.
    # This is the equivalent surface current in the shield that would cause its
    # scalar magnetic potential to be constant
    shield.s = StreamFunction(np.linalg.solve(shield.U_coupling(points), coil.U_coupling(points) @ coil.s),
                              shield)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing scalar potential coupling matrix, 2773 vertices by 2773 target points... took 11.10 seconds.
    Computing scalar potential coupling matrix, 3184 vertices by 2773 target points... took 12.22 seconds.



Plot the difference in field when taking the shield into account


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    B_target = coil.B_coupling(target_points) @ coil.s

    B_target_w_shield = coil.B_coupling(target_points) @ coil.s + shield.B_coupling(target_points) @ shield.s

    B_quiver = mlab.quiver3d(*target_points.T, *(B_target_w_shield - B_target).T, colormap='viridis', mode='arrow')
    f.scene.isometric_view()
    mlab.colorbar(B_quiver, title='Difference in magnetic field (a.u.)')




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_005.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 2773 vertices by 672 target points... took 0.66 seconds.
    This object has no scalar data



Let's redesign the coil taking the shield into account prospectively


.. code-block:: default


    shield.coupling = np.linalg.solve(shield.U_coupling(points), coil.U_coupling(points))

    secondary_C = shield.B_coupling(target_points) @ shield.coupling

    total_C = coil.B_coupling(target_points) + secondary_C

    target_spec_w_shield = {'coupling':total_C, 'rel_error':0, 'abs_error':target_abs_error, 'target':target_field}


    coil.s2, coil.prob2 = optimize_streamfunctions(coil,
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
    Factor     - setup time             : 1.85              dense det. time        : 0.00            
    Factor     - ML order time          : 0.28              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  84.60 
    1   6.6e+01  5.1e-01  9.1e-01  -2.50e-01  1.083065369e+02   1.077565488e+02   5.1e-01  85.17 
    2   3.2e+01  2.4e-01  3.4e-01  -1.80e-01  4.305860933e+02   4.302841974e+02   2.4e-01  85.72 
    3   2.2e+01  1.7e-01  2.0e-01  1.11e+00   6.051999340e+02   6.050078995e+02   1.7e-01  86.23 
    4   8.2e+00  6.4e-02  4.9e-02  9.17e-01   8.680965416e+02   8.680348149e+02   6.4e-02  86.84 
    5   5.7e+00  4.4e-02  3.1e-02  8.95e-01   9.166094221e+02   9.165751401e+02   4.4e-02  87.35 
    6   4.5e-01  3.5e-03  6.7e-04  9.00e-01   1.088049692e+03   1.088046363e+03   3.5e-03  87.96 
    7   2.6e-01  2.0e-03  3.0e-04  9.86e-01   1.095071026e+03   1.095069313e+03   2.0e-03  88.49 
    8   1.3e-01  1.0e-03  1.1e-04  9.92e-01   1.100298352e+03   1.100297589e+03   1.0e-03  89.00 
    9   1.4e-02  1.1e-04  4.0e-06  9.96e-01   1.105213834e+03   1.105213803e+03   1.1e-04  89.61 
    10  1.9e-03  1.5e-05  2.1e-07  1.00e+00   1.105801411e+03   1.105801407e+03   1.5e-05  90.16 
    11  2.3e-04  1.8e-06  7.1e-09  1.00e+00   1.105888972e+03   1.105888972e+03   1.8e-06  90.77 
    12  4.1e-05  3.2e-07  4.6e-10  1.00e+00   1.105898566e+03   1.105898568e+03   3.2e-07  91.33 
    13  2.8e-05  2.2e-07  1.5e-10  1.00e+00   1.105899279e+03   1.105899280e+03   2.2e-07  92.38 
    14  5.5e-06  4.3e-08  1.2e-10  1.00e+00   1.105900448e+03   1.105900449e+03   4.3e-08  92.88 
    15  6.1e-06  2.7e-08  7.2e-11  1.00e+00   1.105900556e+03   1.105900556e+03   2.7e-08  93.84 
    16  6.7e-06  1.6e-08  7.4e-11  1.00e+00   1.105900627e+03   1.105900627e+03   1.6e-08  94.78 
    17  3.0e-06  2.8e-09  2.5e-11  1.00e+00   1.105900716e+03   1.105900715e+03   2.8e-09  95.70 
    Optimizer terminated. Time: 96.17   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.1059007158e+03    nrm: 2e+03    Viol.  con: 5e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.1059007154e+03    nrm: 1e+04    Viol.  con: 6e-08    var: 7e-10    cones: 0e+00  



Plot the newly designed coil windings and field at the target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.s2.vert, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target2 = total_C @ coil.s2
    mlab.quiver3d(*target_points.T, *B_target2.T, mode='arrow', scale_factor=0.75)




    f.scene.isometric_view()
    f.scene.camera.zoom(0.95)





.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_006.png
    :class: sphx-glr-single-img




Plot difference in field


.. code-block:: default



    import seaborn as sns
    import matplotlib.pyplot as plt



    fig, axes = plt.subplots(1, 3, figsize=(12, 3))

    axnames = ['X', 'Y', 'Z']

    #fig.suptitle('Component-wise effect of magnetic shield on target field amplitude distribution')
    for ax_idx, ax in enumerate(axes):

        sns.kdeplot(B_target[:, ax_idx], label='Coil without shield', ax=ax, shade=True, legend=False)
        sns.kdeplot(B_target_w_shield[:, ax_idx], label='Coil with shield', ax=ax, shade=True, legend=False)
        sns.kdeplot(B_target2[:, ax_idx], label='Coil designed with shield', ax=ax, shade=True, legend=False)
    #    ax.set_title(axnames[ax_idx])
        ax.get_yaxis().set_visible(False)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        ax.set_xlabel('Magnetic field on %s-axis'%axnames[ax_idx])

        if ax_idx == 0:
            ax.legend()

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])






.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_007.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 6 minutes  25.552 seconds)

**Estimated memory usage:**  4533 MB


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
