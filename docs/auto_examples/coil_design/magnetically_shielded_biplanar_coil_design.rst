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
    from bfieldtools.magnetic_field_mesh import compute_C, compute_U
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
    shieldmesh = trimesh.load('/l/bfieldtools/bfieldtools/example_meshes/closed_cylinder.stl')
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




Compute C matrices that are used to compute the generated magnetic field


.. code-block:: default


    coil.C = compute_C(coil.mesh, target_points)
    shield.C = compute_C(shield.mesh, target_points)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 3184 vertices by 672 target points... took 0.93 seconds.
    Computing C matrix, 962 vertices by 672 target points... took 0.26 seconds.



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

    target_spec = {'C':coil.C, 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target_field':target_field}

    import mosek

    coil.I, coil.prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 2 chunks since 9 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 68.90 seconds.


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
    Factor     - setup time             : 1.97              dense det. time        : 0.00            
    Factor     - ML order time          : 0.32              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  77.63 
    1   4.0e+01  6.2e-01  2.2e-01  1.05e+00   4.677130957e+01   4.602040536e+01   6.2e-01  78.31 
    2   8.8e+00  1.4e-01  2.0e-02  1.19e+00   7.851531865e+01   7.839683819e+01   1.4e-01  78.94 
    3   4.2e+00  6.6e-02  9.0e-03  1.30e+00   8.060276587e+01   8.055197009e+01   6.6e-02  79.54 
    4   3.2e+00  5.0e-02  7.0e-03  1.01e+00   8.118900965e+01   8.114766773e+01   5.0e-02  80.19 
    5   9.6e-02  1.5e-03  3.8e-05  1.09e+00   8.540388976e+01   8.540282440e+01   1.5e-03  80.87 
    6   1.0e-02  1.6e-04  1.7e-06  1.01e+00   8.545287301e+01   8.545276641e+01   1.6e-04  81.62 
    7   1.1e-03  1.7e-05  5.6e-08  1.00e+00   8.547050333e+01   8.547049238e+01   1.7e-05  82.65 
    8   2.2e-04  3.5e-06  5.3e-09  9.99e-01   8.547222412e+01   8.547222185e+01   3.5e-06  83.78 
    9   1.1e-04  1.7e-06  1.9e-09  1.00e+00   8.547245153e+01   8.547245040e+01   1.7e-06  84.69 
    10  4.2e-06  6.6e-08  1.4e-11  1.00e+00   8.547266623e+01   8.547266619e+01   6.6e-08  86.07 
    11  8.4e-07  1.6e-08  1.2e-11  1.00e+00   8.547267483e+01   8.547267416e+01   3.3e-10  87.22 
    Optimizer terminated. Time: 87.84   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 8.5472674826e+01    nrm: 2e+02    Viol.  con: 1e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 8.5472674165e+01    nrm: 1e+03    Viol.  con: 1e-08    var: 1e-09    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    mlab.quiver3d(*target_points.T, *B_target.T)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_002.png
    :class: sphx-glr-single-img




Now, let's compute the effect of the shield on the field produced by the coil


.. code-block:: default


    # Calculate primary potential matrix at the shield surface
    P_prim = compute_U(coil.mesh, shield.mesh.vertices)

    # Calculate linear collocation BEM matrix
    P_bem = compute_U(shield.mesh, shield.mesh.vertices)

    # Recalculate diag elements according to de Munck paper
    for diag_index in range(P_bem.shape[0]):
        P_bem[diag_index, diag_index] = 0
        P_bem[diag_index, diag_index] = -P_bem[diag_index, :].sum()

    # Matrix misses one rank, make it invertible
    # by rank-one update (sets potential of constant dipole layer)
    P_bem += np.ones(P_bem.shape)/P_bem.shape[0]


    # Solve equivalent stream function for the perfect linear mu-metal layer.
    # This is the equivalent surface current in the shield that would cause its
    # scalar magnetic potential to be constant
    shield.I =  np.linalg.solve(P_bem, P_prim @ coil.I)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing U matrix, 3184 vertices by 962 target points... took 13.88 seconds.
    Computing U matrix, 962 vertices by 962 target points... took 4.01 seconds.



Plot the difference in field when taking the shield into account


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    B_target_w_shield = coil.C.transpose([0, 2, 1]) @ coil.I + shield.C.transpose([0, 2, 1]) @ shield.I

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

    This object has no scalar data



Let's redesign the coil taking the shield into account prospectively


.. code-block:: default


    shield.coupling = np.linalg.solve(P_bem, P_prim)

    secondary_C = (shield.C.transpose((0,2,1)) @ shield.coupling).transpose((0,2,1))

    total_C = coil.C + secondary_C

    target_spec_w_shield = {'C':total_C, 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target_field':target_field}


    coil.I2, coil.prob2 = optimize_streamfunctions(coil,
                                       [target_spec_w_shield],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none



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
    Factor     - setup time             : 1.88              dense det. time        : 0.00            
    Factor     - ML order time          : 0.28              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 4.93e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  68.36 
    1   3.7e+01  5.8e-01  3.6e-01  1.09e+00   4.344673635e+01   4.275006379e+01   5.8e-01  68.96 
    2   5.9e+00  9.1e-02  1.5e-02  1.23e+00   5.743682235e+01   5.735225303e+01   9.1e-02  69.60 
    3   2.3e+00  3.7e-02  5.0e-03  1.23e+00   5.886873114e+01   5.883918399e+01   3.7e-02  70.20 
    4   7.7e-01  1.2e-02  1.0e-03  1.09e+00   5.952390006e+01   5.951485524e+01   1.2e-02  70.77 
    5   6.1e-01  9.4e-03  7.3e-04  1.02e+00   5.962053624e+01   5.961346752e+01   9.4e-03  71.33 
    6   9.1e-02  1.4e-03  4.8e-05  1.02e+00   5.986565873e+01   5.986466682e+01   1.4e-03  72.09 
    7   7.4e-03  1.1e-04  1.1e-06  9.99e-01   5.993844588e+01   5.993836570e+01   1.1e-04  72.84 
    8   5.5e-03  8.6e-05  7.1e-07  1.00e+00   5.994005206e+01   5.993999214e+01   8.6e-05  73.43 
    9   4.3e-03  6.7e-05  5.0e-07  9.99e-01   5.994130700e+01   5.994125973e+01   6.7e-05  74.00 
    10  1.4e-03  2.2e-05  9.0e-08  1.00e+00   5.994398039e+01   5.994396533e+01   2.2e-05  74.56 
    11  1.5e-04  2.3e-06  3.1e-09  1.00e+00   5.994529199e+01   5.994529038e+01   2.3e-06  75.17 
    12  3.7e-06  5.7e-08  8.5e-12  1.00e+00   5.994544991e+01   5.994544988e+01   5.7e-08  75.81 
    13  1.8e-06  2.9e-08  1.1e-12  1.00e+00   5.994545192e+01   5.994545193e+01   2.9e-08  76.91 
    14  9.2e-07  1.4e-08  8.2e-13  1.00e+00   5.994545292e+01   5.994545292e+01   1.4e-08  77.79 
    Optimizer terminated. Time: 78.28   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 5.9945452924e+01    nrm: 1e+02    Viol.  con: 6e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 5.9945452920e+01    nrm: 9e+02    Viol.  con: 1e-06    var: 4e-10    cones: 0e+00  



Plot the newly designed coil windings and field at the target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.I2, N_contours=10)
    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target2 = total_C.transpose([0, 2, 1]) @ coil.I2
    mlab.quiver3d(*target_points.T, *B_target2.T)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_005.png
    :class: sphx-glr-single-img




Plot the difference in stream functions


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_data_on_vertices(coil.mesh, np.nan_to_num(100 * (coil.I-coil.I2)/coil.I), figure=f, colorbar=True)

    mlab.colorbar(title='Relative error (%)')





.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_006.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /l/bfieldtools/examples/coil_design/magnetically_shielded_biplanar_coil_design.py:244: RuntimeWarning: invalid value encountered in true_divide
      plot_data_on_vertices(coil.mesh, np.nan_to_num(100 * (coil.I-coil.I2)/coil.I), figure=f, colorbar=True)



Finally, plot the field lines when the shield is included into the model


.. code-block:: default


    extent = 8
    N = 20
    X, Y, Z = np.meshgrid(np.linspace(-extent, extent, N)+7.5, np.linspace(-extent, extent, N), np.linspace(-extent, extent, N))

    r = np.array([X.flatten(), Y.flatten(), Z.flatten()]).T

    r = r[shield.mesh.contains(r)]


    coil.C_cyl = compute_C(coil.mesh, r)
    shield.C_cyl = compute_C(shield.mesh, r)

    secondary_C_cyl = (shield.C_cyl.transpose((0,2,1)) @ shield.coupling).transpose((0,2,1))

    total_C_cyl = coil.C_cyl + secondary_C_cyl


    Bfield = total_C_cyl.transpose([0, 2, 1]) @ coil.I2

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    quiv = mlab.quiver3d(*r.T, *Bfield.T)



    plot_3d_current_loops(loops, colors='auto', figure=f)

    shield.plot_mesh(representation='surface', opacity=0.1, cull_front=True)



.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_007.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 3184 vertices by 4856 target points... took 5.18 seconds.
    Computing C matrix, 962 vertices by 4856 target points... took 1.56 seconds.




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 7 minutes  9.468 seconds)

**Estimated memory usage:**  7947 MB


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
