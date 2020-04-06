.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_coil_design_magnetically_shielded_biplanar_coil_design.py>`     to download the full example code
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


    from bfieldtools.conductor import Conductor, StreamFunction
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

    Computing magnetic field coupling matrix, 3184 vertices by 672 target points... took 0.96 seconds.
    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 80 chunks (9723 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 56.47 seconds.
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
    Factor     - setup time             : 1.63              dense det. time        : 0.00            
    Factor     - ML order time          : 0.19              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  104.97
    1   6.2e+01  4.9e-01  8.3e-01  -1.88e-01  1.238067944e+02   1.232654642e+02   4.9e-01  106.36
    2   2.8e+01  2.2e-01  2.7e-01  -7.61e-02  4.720693519e+02   4.717781366e+02   2.2e-01  107.66
    3   6.6e+00  5.2e-02  3.4e-02  1.16e+00   7.660921118e+02   7.660494017e+02   5.2e-02  109.08
    4   1.2e+00  9.7e-03  2.7e-03  1.01e+00   8.627082696e+02   8.626988383e+02   9.7e-03  110.58
    5   1.7e-01  1.3e-03  1.4e-04  9.78e-01   8.843878725e+02   8.843868307e+02   1.3e-03  111.91
    6   2.1e-02  1.6e-04  6.4e-06  1.00e+00   8.876841755e+02   8.876840522e+02   1.6e-04  113.25
    7   2.6e-03  2.0e-05  2.8e-07  1.00e+00   8.881627682e+02   8.881627508e+02   2.0e-05  114.55
    8   3.6e-04  2.8e-06  1.4e-08  1.00e+00   8.882189076e+02   8.882189053e+02   2.8e-06  115.78
    9   1.7e-07  1.0e-07  1.8e-12  1.00e+00   8.882280594e+02   8.882280594e+02   1.0e-09  117.17
    10  1.2e-05  5.1e-08  1.3e-12  1.00e+00   8.882280608e+02   8.882280608e+02   5.0e-10  118.94
    11  2.1e-05  2.6e-08  6.8e-12  1.00e+00   8.882280616e+02   8.882280620e+02   2.5e-10  120.38
    12  5.0e-05  1.3e-08  2.1e-12  1.00e+00   8.882280619e+02   8.882280617e+02   1.3e-10  121.83
    Optimizer terminated. Time: 122.34  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 8.8822806192e+02    nrm: 2e+03    Viol.  con: 2e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 8.8822806170e+02    nrm: 6e+03    Viol.  con: 6e-09    var: 5e-09    cones: 2e-13  




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

    Computing scalar potential coupling matrix, 2773 vertices by 2773 target points... took 14.56 seconds.
    Computing scalar potential coupling matrix, 3184 vertices by 2773 target points... took 16.14 seconds.




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

    Computing magnetic field coupling matrix, 2773 vertices by 672 target points... took 0.90 seconds.
    This object has no scalar data

    <mayavi.core.lut_manager.LUTManager object at 0x000001928787FC50>



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
    Factor     - setup time             : 1.50              dense det. time        : 0.00            
    Factor     - ML order time          : 0.19              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  105.55
    1   6.6e+01  5.1e-01  9.2e-01  -2.50e-01  1.083155072e+02   1.077649901e+02   5.1e-01  106.89
    2   3.2e+01  2.5e-01  3.4e-01  -1.80e-01  4.307431787e+02   4.304407604e+02   2.5e-01  108.13
    3   2.2e+01  1.7e-01  2.0e-01  1.11e+00   6.039188441e+02   6.037257037e+02   1.7e-01  109.42
    4   8.2e+00  6.4e-02  4.9e-02  9.17e-01   8.682557774e+02   8.681939527e+02   6.4e-02  110.75
    5   5.7e+00  4.4e-02  3.1e-02  8.95e-01   9.169496595e+02   9.169153545e+02   4.4e-02  112.00
    6   4.3e-01  3.3e-03  6.3e-04  8.99e-01   1.089000379e+03   1.088997032e+03   3.3e-03  113.31
    7   2.5e-01  1.9e-03  2.8e-04  9.87e-01   1.095898711e+03   1.095897014e+03   1.9e-03  114.61
    8   1.3e-01  9.7e-04  1.0e-04  9.93e-01   1.100927626e+03   1.100926886e+03   9.7e-04  115.91
    9   1.4e-02  1.1e-04  4.3e-06  9.96e-01   1.105489681e+03   1.105489650e+03   1.1e-04  117.25
    10  1.7e-03  1.3e-05  1.7e-07  1.00e+00   1.106121892e+03   1.106121889e+03   1.3e-05  118.52
    11  1.9e-05  1.8e-07  5.7e-11  1.00e+00   1.106207121e+03   1.106207121e+03   1.5e-07  119.97
    12  1.4e-05  1.3e-07  6.0e-11  1.00e+00   1.106207405e+03   1.106207406e+03   1.1e-07  121.58
    13  4.5e-06  4.1e-08  5.1e-13  1.00e+00   1.106207874e+03   1.106207874e+03   3.5e-08  122.83
    14  3.0e-06  2.8e-08  4.7e-12  1.00e+00   1.106207952e+03   1.106207952e+03   2.4e-08  124.42
    15  9.5e-06  6.3e-09  7.2e-12  1.00e+00   1.106208075e+03   1.106208075e+03   5.2e-09  125.94
    Optimizer terminated. Time: 126.44  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.1062080747e+03    nrm: 2e+03    Viol.  con: 1e-08    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.1062080746e+03    nrm: 1e+04    Viol.  con: 1e-07    var: 4e-09    cones: 0e+00  




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

   **Total running time of the script:** ( 7 minutes  14.488 seconds)


.. _sphx_glr_download_auto_examples_coil_design_magnetically_shielded_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: magnetically_shielded_biplanar_coil_design.py <magnetically_shielded_biplanar_coil_design.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: magnetically_shielded_biplanar_coil_design.ipynb <magnetically_shielded_biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
