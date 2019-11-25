.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_biplanar_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_biplanar_coil_design.py:


Biplanar coil design
====================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes.



.. code-block:: default


    import numpy as np
    import matplotlib.pyplot as plt
    from mayavi import mlab
    import trimesh

    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

    planemesh.apply_scale(scaling_factor*1.6)

    #Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 5, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                             planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                         planemesh.faces, process=False)

    joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = MeshWrapper(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)







Set up target and stray field points


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 0]) * scaling_factor

    sidelength = 2 * scaling_factor
    n = 8
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



    #    #Here, the stray field points are on a spherical surface
    stray_radius = 20 * scaling_factor
    #    stray_length = 20 * scaling_factor
    #
    #    stray_points = cylinder_points(radius=stray_radius,
    #                                   length = stray_length,
    #                                   nlength = 5,
    #                                   nalpha = 30,
    #                                   orientation=np.array([1, 0, 0]))
    #
    stray_points_mesh = trimesh.creation.icosphere(subdivisions=3, radius=stray_radius)
    stray_points = stray_points_mesh.vertices + center

    n_stray_points = len(stray_points)








Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function

    target_field = np.zeros(target_points.shape)
    target_field[:, 0] = target_field[:, 0] + 1

    target_rel_error = np.zeros_like(target_field)
    target_rel_error[:, 0] += 0.01

    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 0] += 0.001
    target_abs_error[:, 1:3] += 0.005

    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}
    stray_spec = {'coupling':coil.B_coupling(stray_points), 'abs_error':0.01, 'rel_error':0, 'target':np.zeros((n_stray_points, 3))}

    bfield_specification = [target_spec, stray_spec]





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 3184 vertices by 160 target points... took 0.31 seconds.
    Computing C matrix, 3184 vertices by 642 target points... took 0.93 seconds.



Run QP solver


.. code-block:: default

    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       [target_spec, stray_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 2 chunks since 7 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 73.82 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 7710            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 7710            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 7710              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 2.01              dense det. time        : 0.00            
    Factor     - ML order time          : 0.28              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 5.26e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.4e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  73.45 
    1   3.0e+01  4.7e-01  1.2e+00  -8.53e-01  6.713334268e+00   6.642742483e+00   4.7e-01  74.12 
    2   2.2e+01  3.4e-01  9.5e-01  -6.35e-01  2.340534155e+01   2.373451278e+01   3.4e-01  74.67 
    3   1.3e+01  2.0e-01  6.1e-01  -4.98e-01  1.132349243e+02   1.140777096e+02   2.0e-01  75.23 
    4   8.4e+00  1.3e-01  3.9e-01  -2.63e-01  1.466684480e+02   1.477520639e+02   1.3e-01  75.80 
    5   3.1e+00  4.8e-02  1.2e-01  -1.99e-02  3.686600597e+02   3.697121573e+02   4.8e-02  76.36 
    6   1.3e+00  2.0e-02  3.6e-02  4.97e-01   4.480844978e+02   4.486499110e+02   2.0e-02  76.97 
    7   1.0e+00  1.6e-02  2.8e-02  6.73e-01   4.573715522e+02   4.579365541e+02   1.6e-02  77.56 
    8   7.1e-01  1.1e-02  1.9e-02  5.06e-01   5.658891744e+02   5.664395502e+02   1.1e-02  78.13 
    9   5.0e-01  7.8e-03  1.2e-02  4.88e-01   8.450951961e+02   8.455168378e+02   7.8e-03  78.72 
    10  4.0e-01  6.2e-03  9.4e-03  5.32e-01   9.163853989e+02   9.168145036e+02   6.2e-03  79.30 
    11  1.1e-01  1.7e-03  1.6e-03  6.68e-01   1.486504636e+03   1.486670483e+03   1.7e-03  80.08 
    12  5.4e-02  8.4e-04  6.2e-04  7.90e-01   1.741137001e+03   1.741245548e+03   8.4e-04  80.64 
    13  7.6e-03  1.2e-04  3.8e-05  9.06e-01   1.967976255e+03   1.967998803e+03   1.2e-04  81.25 
    14  8.8e-04  1.4e-05  1.5e-06  9.79e-01   2.018240273e+03   2.018243006e+03   1.4e-05  81.84 
    15  2.4e-06  1.6e-07  2.1e-10  9.97e-01   2.025562893e+03   2.025562900e+03   3.8e-08  82.67 
    16  1.1e-06  1.4e-07  1.2e-10  1.00e+00   2.025573181e+03   2.025573186e+03   1.8e-08  83.79 
    17  2.6e-07  3.2e-08  4.2e-11  1.00e+00   2.025580480e+03   2.025580481e+03   4.0e-09  84.88 
    18  2.5e-07  3.1e-08  3.8e-11  1.00e+00   2.025580537e+03   2.025580538e+03   3.9e-09  86.09 
    19  2.5e-07  3.1e-08  1.5e-11  1.00e+00   2.025580566e+03   2.025580566e+03   3.9e-09  87.29 
    20  2.4e-07  3.1e-08  1.9e-11  1.00e+00   2.025580594e+03   2.025580595e+03   3.8e-09  88.47 
    21  4.5e-07  3.9e-09  1.0e-11  1.00e+00   2.025582365e+03   2.025582366e+03   4.7e-10  89.52 
    22  4.5e-07  3.9e-09  1.0e-11  1.00e+00   2.025582365e+03   2.025582366e+03   4.7e-10  90.88 
    23  4.5e-07  3.9e-09  1.0e-11  1.00e+00   2.025582365e+03   2.025582366e+03   4.7e-10  92.20 
    Optimizer terminated. Time: 94.09   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 2.0255823651e+03    nrm: 4e+03    Viol.  con: 7e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 2.0255823657e+03    nrm: 9e+04    Viol.  con: 2e-05    var: 1e-08    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)








.. image:: /auto_examples/coil_design/images/sphx_glr_biplanar_coil_design_001.png
    :class: sphx-glr-single-img




Plot cross-section of magentic field and magnetic potential of the discretized loops


.. code-block:: default


    from bfieldtools.line_magnetics import magnetic_field, scalar_potential

    x = y = np.linspace(-12, 12, 250)
    X,Y = np.meshgrid(x, y, indexing='ij')
    points = np.zeros((X.flatten().shape[0], 3))
    points[:, 0] = X.flatten()
    points[:, 1] = Y.flatten()

    B = np.zeros_like(points)
    U = np.zeros((points.shape[0],))
    for loop_idx in range(len(loops)):
        B += magnetic_field(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)
        U += scalar_potential(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)

    B = B.T[:2].reshape(2, x.shape[0], y.shape[0])
    lw = np.sqrt(B[0]**2 + B[1]**2)
    lw = 2*lw/np.max(lw)

    U = U.reshape(x.shape[0], y.shape[0])

    plt.figure()

    plt.pcolormesh(X, Y, U.T, cmap='seismic', shading='gouraud')
    #plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
    #           extent=(x.min(), x.max(), y.min(), y.max()))

    seed_points=points[:,:2]*0.3

    plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k', integration_direction='both',
                   start_points=seed_points)
    plt.axis('equal')
    plt.axis('off')
    for loop in loops:
        plt.plot(loop[:,1], loop[:,0], 'k', linewidth=4, alpha=0.1)

    plt.tight_layout()



.. image:: /auto_examples/coil_design/images/sphx_glr_biplanar_coil_design_002.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 7 minutes  35.863 seconds)

**Estimated memory usage:**  8011 MB


.. _sphx_glr_download_auto_examples_coil_design_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: biplanar_coil_design.py <biplanar_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: biplanar_coil_design.ipynb <biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
