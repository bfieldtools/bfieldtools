.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_spherical_harmonics_biplanar_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_spherical_harmonics_biplanar_coil_design.py:


High-order spherical harmonic biplanar coil design
==================================================

Example showing a basic biplanar coil producing a high-order spherical harmonic field
in a target region between the two coil planes.



.. code-block:: default


    import numpy as np
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











Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function


    from bfieldtools.sphtools import sphbasis


    sph = sphbasis(50)

    #plotsph.plotYlms(sph, 3)

    lmax = 4
    alm = np.zeros((lmax*(lmax+2),))
    blm = np.zeros((lmax*(lmax+2),))

    #
    #alm[22]+=1
    blm[22]+=1

    sphfield = sph.field(target_points, alm, blm, lmax)

    target_field = sphfield/np.max(sphfield[:, 0])

    #target_field[:, 2] = 0


    coil.plot_mesh()
    mlab.quiver3d(*target_points.T, *sphfield.T)



    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':0, 'abs_error':0.1, 'target':target_field}
    stray_spec = {'coupling':coil.B_coupling(stray_points), 'abs_error':0.01, 'rel_error':0, 'target':np.zeros((n_stray_points, 3))}

    bfield_specification = [target_spec, stray_spec]




.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_biplanar_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.29 seconds.
    Computing magnetic field coupling matrix, 3184 vertices by 642 target points... took 0.88 seconds.



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

    Computing inductance matrix in 2 chunks since 8 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 72.04 seconds.
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
    Factor     - setup time             : 2.00              dense det. time        : 0.00            
    Factor     - ML order time          : 0.28              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 5.26e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   3.3e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  80.40 
    1   2.6e+02  8.0e-01  1.8e+00  -9.69e-01  1.155643995e+00   3.990860718e-01   8.0e-01  81.01 
    2   2.0e+02  6.2e-01  1.6e+00  -9.57e-01  1.172923928e+01   1.130691247e+01   6.2e-01  81.58 
    3   1.8e+02  5.4e-01  1.4e+00  -9.42e-01  5.150142546e+01   5.130930199e+01   5.4e-01  82.13 
    4   1.5e+02  4.4e-01  1.3e+00  -9.36e-01  5.688278369e+01   5.704030827e+01   4.4e-01  82.69 
    5   1.2e+02  3.6e-01  1.1e+00  -9.26e-01  3.706864434e+02   3.713022974e+02   3.6e-01  83.25 
    6   2.7e+01  8.4e-02  4.6e-01  -9.00e-01  4.514670216e+03   4.520531449e+03   8.4e-02  83.95 
    7   6.7e+00  2.0e-02  1.8e-01  -7.82e-01  2.322274901e+04   2.324145982e+04   2.0e-02  84.54 
    8   1.4e+00  4.2e-03  4.2e-02  -4.03e-01  7.606520131e+04   7.608951558e+04   4.2e-03  85.09 
    9   5.6e-01  1.7e-03  1.2e-02  3.48e-01   1.008058436e+05   1.008181630e+05   1.7e-03  85.66 
    10  1.8e-02  5.6e-05  9.0e-05  7.19e-01   1.223248483e+05   1.223254952e+05   5.6e-05  86.47 
    11  2.9e-03  8.9e-06  5.2e-06  9.88e-01   1.230883141e+05   1.230883996e+05   8.9e-06  87.06 
    12  1.4e-05  4.1e-08  4.3e-08  9.98e-01   1.232536590e+05   1.232536596e+05   4.2e-08  87.78 
    13  6.8e-06  2.1e-08  1.9e-08  1.00e+00   1.232540549e+05   1.232540553e+05   2.1e-08  88.99 
    14  8.0e-06  1.9e-08  1.6e-08  1.00e+00   1.232540798e+05   1.232540799e+05   2.0e-08  90.07 
    15  5.6e-06  1.7e-08  2.7e-09  1.00e+00   1.232541264e+05   1.232541266e+05   1.7e-08  91.10 
    16  5.6e-06  1.7e-08  2.7e-09  1.00e+00   1.232541264e+05   1.232541266e+05   1.7e-08  92.36 
    17  5.6e-06  1.7e-08  2.7e-09  1.00e+00   1.232541264e+05   1.232541266e+05   1.7e-08  93.63 
    18  2.8e-06  8.5e-09  1.0e-09  1.00e+00   1.232542896e+05   1.232542897e+05   8.6e-09  94.71 
    19  2.8e-06  8.5e-09  1.0e-09  1.00e+00   1.232542896e+05   1.232542897e+05   8.6e-09  95.99 
    20  3.9e-06  8.4e-09  5.7e-10  1.00e+00   1.232542921e+05   1.232542922e+05   8.4e-09  97.02 
    21  3.9e-06  8.4e-09  5.7e-10  1.00e+00   1.232542921e+05   1.232542922e+05   8.4e-09  98.31 
    22  3.9e-06  8.4e-09  5.7e-10  1.00e+00   1.232542921e+05   1.232542922e+05   8.4e-09  99.52 
    Optimizer terminated. Time: 101.26  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.2325429211e+05    nrm: 2e+05    Viol.  con: 9e-07    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.2325429219e+05    nrm: 5e+05    Viol.  con: 2e-02    var: 2e-07    cones: 0e+00  



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


.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_biplanar_coil_design_002.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 3 minutes  48.448 seconds)

**Estimated memory usage:**  7986 MB


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
