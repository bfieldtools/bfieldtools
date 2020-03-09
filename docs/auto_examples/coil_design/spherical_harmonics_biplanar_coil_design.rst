.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_spherical_harmonics_biplanar_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_spherical_harmonics_biplanar_coil_design.py:


High-order spherical harmonic biplanar coil design
==================================================

Example showing a basic biplanar coil producing a high-order spherical harmonic field
in a specific target region between the two coil planes.



.. code-block:: default


    import numpy as np
    from mayavi import mlab
    import trimesh

    from bfieldtools.mesh_class import Conductor
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
    coil = Conductor(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)







Set up target and stray field points


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 0]) * scaling_factor

    sidelength = 2 * scaling_factor
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


    from bfieldtools import sphtools


    lmax = 4
    alm = np.zeros((lmax*(lmax+2),))
    blm = np.zeros((lmax*(lmax+2),))

    blm[22]+=1

    sphfield = sphtools.field(target_points, alm, blm, lmax)

    target_field = sphfield/np.max(sphfield[:, 0])



    coil.plot_mesh(opacity=0.2)
    mlab.quiver3d(*target_points.T, *sphfield.T)



    target_spec = {'coupling':coil.B_coupling(target_points), 'abs_error':0.1, 'target':target_field}
    stray_spec = {'coupling':coil.B_coupling(stray_points), 'abs_error':0.01, 'target':np.zeros((n_stray_points, 3))}

    bfield_specification = [target_spec, stray_spec]




.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_biplanar_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 672 target points... took 0.76 seconds.
    Computing magnetic field coupling matrix, 3184 vertices by 642 target points... took 0.69 seconds.



Run QP solver


.. code-block:: default

    import mosek

    coil.s, prob = optimize_streamfunctions(coil,
                                       [target_spec, stray_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 34964 MiB required for 3184 by 3184 vertices...
    Computing inductance matrix in 60 chunks (11705 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 42.45 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 10782           
      Cones                  : 1               
      Scalar variables       : 6083            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 10782           
      Cones                  : 1               
      Scalar variables       : 6083            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 10782             conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 2.82              dense det. time        : 0.00            
    Factor     - ML order time          : 0.30              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 6.55e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.6e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  193.66
    1   5.2e+02  8.0e-01  1.8e+00  -9.86e-01  4.923004185e+00   4.170266302e+00   8.0e-01  194.39
    2   3.6e+02  5.6e-01  1.5e+00  -9.81e-01  5.271060710e+01   5.249286736e+01   5.6e-01  195.08
    3   2.7e+02  4.1e-01  1.3e+00  -9.69e-01  2.410919996e+02   2.414514590e+02   4.1e-01  195.76
    4   2.2e+02  3.3e-01  1.1e+00  -9.63e-01  2.680807755e+02   2.690147129e+02   3.3e-01  196.44
    5   1.5e+02  2.3e-01  9.2e-01  -9.55e-01  1.824694249e+03   1.826794837e+03   2.3e-01  197.13
    6   2.1e+01  3.2e-02  2.9e-01  -9.26e-01  2.544651385e+04   2.546414809e+04   3.2e-02  198.11
    7   9.5e+00  1.4e-02  1.7e-01  -7.56e-01  6.646639493e+04   6.649953268e+04   1.4e-02  198.82
    8   7.8e+00  1.2e-02  1.5e-01  -5.96e-01  8.105878983e+04   8.109584779e+04   1.2e-02  199.51
    9   1.0e+00  1.5e-03  2.9e-02  -5.92e-01  4.395893298e+05   4.396771666e+05   1.5e-03  200.48
    10  5.3e-01  8.1e-04  1.2e-02  2.98e-01   6.131845050e+05   6.132416891e+05   8.1e-04  201.16
    11  4.2e-01  6.4e-04  9.8e-03  2.64e-01   6.481199127e+05   6.481784877e+05   6.4e-04  201.85
    12  1.3e-01  2.1e-04  1.9e-03  6.19e-01   8.432815698e+05   8.433029120e+05   2.1e-04  202.53
    13  3.7e-03  5.6e-06  1.1e-05  8.51e-01   9.493346919e+05   9.493355611e+05   5.6e-06  203.38
    14  5.8e-04  2.7e-06  6.7e-07  9.94e-01   9.523270588e+05   9.523271899e+05   8.9e-07  204.11
    15  3.0e-04  1.4e-06  4.1e-07  9.99e-01   9.526052592e+05   9.526053261e+05   4.5e-07  205.73
    16  5.7e-06  2.1e-08  3.0e-08  1.00e+00   9.528934259e+05   9.528934272e+05   6.0e-09  206.94
    17  5.7e-06  2.1e-08  2.9e-08  1.00e+00   9.528934525e+05   9.528934538e+05   5.9e-09  208.35
    18  3.9e-05  1.2e-08  9.3e-09  1.00e+00   9.528951446e+05   9.528951449e+05   3.3e-09  209.55
    19  5.2e-05  1.0e-08  7.9e-09  1.00e+00   9.528953874e+05   9.528953880e+05   3.0e-09  210.82
    20  6.5e-05  9.1e-09  2.5e-09  1.00e+00   9.528956041e+05   9.528956044e+05   2.6e-09  212.10
    21  6.4e-05  9.0e-09  5.7e-10  1.00e+00   9.528956162e+05   9.528956165e+05   2.6e-09  213.49
    22  6.4e-05  9.0e-09  5.7e-10  1.00e+00   9.528956162e+05   9.528956165e+05   2.6e-09  215.07
    23  6.4e-05  9.0e-09  5.7e-10  1.00e+00   9.528956162e+05   9.528956165e+05   2.6e-09  216.67
    Optimizer terminated. Time: 218.97  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 9.5289561616e+05    nrm: 2e+06    Viol.  con: 1e-06    var: 0e+00    cones: 0e+00  
      Dual.    obj: 9.5289561653e+05    nrm: 1e+07    Viol.  con: 3e-02    var: 9e-07    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.s, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.s

    mlab.quiver3d(*target_points.T, *B_target.T)


.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_biplanar_coil_design_002.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 6 minutes  12.030 seconds)

**Estimated memory usage:**  6654 MB


.. _sphx_glr_download_auto_examples_coil_design_spherical_harmonics_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: spherical_harmonics_biplanar_coil_design.py <spherical_harmonics_biplanar_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: spherical_harmonics_biplanar_coil_design.ipynb <spherical_harmonics_biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
