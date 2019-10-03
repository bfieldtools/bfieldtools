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
    from bfieldtools.magnetic_field_mesh import compute_C
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

    planemesh.apply_scale(scaling_factor)

    #Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 4, 0]) * scaling_factor

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









Compute C matrices that are used to compute the generated magnetic field


.. code-block:: default


    coil.C = compute_C(coil.mesh, target_points)
    coil.strayC = compute_C(coil.mesh, stray_points)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 3184 vertices by 160 target points... took 0.30 seconds.
    Computing C matrix, 3184 vertices by 642 target points... took 0.84 seconds.



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

    target_spec = {'C':coil.C, 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target_field':target_field}
    stray_spec = {'C':coil.strayC, 'abs_error':0.01, 'rel_error':0, 'target_field':np.zeros((n_stray_points, 3))}

    bfield_specification = [target_spec, stray_spec]







Run QP solver


.. code-block:: default

    import mosek

    coil.I, prob = optimize_streamfunctions(coil,
                                       [target_spec, stray_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 2 chunks since 8 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 65.46 seconds.


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
    Factor     - ML order time          : 0.29              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 5.26e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   3.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  75.60 
    1   1.8e+02  5.5e-01  1.5e+00  -9.69e-01  3.901237432e-01   1.731532757e-01   5.5e-01  76.35 
    2   9.6e+01  2.9e-01  1.0e+00  -9.27e-01  4.185471306e+00   5.344659909e+00   2.9e-01  76.89 
    3   8.0e+01  2.4e-01  9.1e-01  -8.52e-01  7.541655539e+00   9.185280969e+00   2.4e-01  77.43 
    4   2.7e+01  8.2e-02  4.6e-01  -8.22e-01  6.067138585e+01   6.676969596e+01   8.2e-02  78.15 
    5   1.8e+01  5.5e-02  3.2e-01  -5.35e-01  1.129291606e+02   1.202495325e+02   5.5e-02  78.73 
    6   3.2e+00  9.7e-03  6.3e-02  -3.68e-01  4.183598608e+02   4.284717884e+02   9.7e-03  79.34 
    7   9.6e-01  2.9e-03  1.0e-02  5.15e-01   5.931225935e+02   5.961146234e+02   2.9e-03  79.87 
    8   4.3e-01  1.3e-03  3.1e-03  9.19e-01   4.222103767e+02   4.235030322e+02   1.3e-03  80.44 
    9   3.3e-01  1.0e-03  2.1e-03  9.88e-01   4.916049319e+02   4.926280323e+02   1.0e-03  80.98 
    10  8.3e-03  2.5e-05  7.2e-06  9.94e-01   4.664898365e+02   4.665081857e+02   2.5e-05  81.73 
    11  8.3e-04  2.5e-06  2.8e-07  1.00e+00   4.673976385e+02   4.674004349e+02   2.5e-06  82.48 
    12  5.1e-04  1.5e-06  1.3e-07  1.00e+00   4.674257060e+02   4.674274237e+02   1.5e-06  83.19 
    13  1.5e-04  4.6e-07  2.1e-08  1.00e+00   4.675038693e+02   4.675043792e+02   4.6e-07  83.81 
    14  1.9e-05  5.8e-08  9.7e-10  1.00e+00   4.675437418e+02   4.675438062e+02   5.8e-08  84.47 
    15  7.3e-06  2.2e-08  2.3e-10  1.00e+00   4.675482297e+02   4.675482544e+02   2.2e-08  85.03 
    16  1.0e-06  5.7e-09  1.4e-11  1.00e+00   4.675502446e+02   4.675502480e+02   3.1e-09  85.59 
    17  2.1e-08  1.2e-10  5.4e-14  1.00e+00   4.675506350e+02   4.675506350e+02   6.5e-11  86.72 
    Optimizer terminated. Time: 87.23   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 4.6755063496e+02    nrm: 9e+02    Viol.  con: 2e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 4.6755063503e+02    nrm: 4e+04    Viol.  con: 7e-05    var: 1e-09    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    mlab.quiver3d(*target_points.T, *B_target.T)


.. image:: /auto_examples/coil_design/images/sphx_glr_biplanar_coil_design_001.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 3 minutes  21.344 seconds)

**Estimated memory usage:**  7944 MB


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
