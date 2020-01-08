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

    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.27 seconds.
    Computing magnetic field coupling matrix, 3184 vertices by 642 target points... took 0.72 seconds.



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

    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 405514 MiB required for 3184 times 3184 vertices...
    Computing inductance matrix in 43 chunks since 9459 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 87.61 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


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
    Factor     - setup time             : 1.83              dense det. time        : 0.00            
    Factor     - ML order time          : 0.16              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 5.26e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.6e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  100.16
    1   9.1e+01  5.6e-01  1.4e+00  -9.59e-01  7.080359218e+00   6.819759266e+00   5.6e-01  100.84
    2   8.0e+01  4.9e-01  1.3e+00  -9.11e-01  3.678727056e+01   3.674375251e+01   4.9e-01  101.38
    3   6.9e+01  4.2e-01  1.2e+00  -9.12e-01  4.814197203e+01   4.836534158e+01   4.2e-01  101.92
    4   5.9e+01  3.6e-01  1.1e+00  -9.02e-01  4.011336960e+02   4.016798419e+02   3.6e-01  102.46
    5   4.3e+01  2.6e-01  9.4e-01  -8.74e-01  4.582647612e+02   4.596444320e+02   2.6e-01  103.00
    6   7.4e+00  4.5e-02  3.3e-01  -8.98e-01  9.180843002e+03   9.192794388e+03   4.5e-02  103.77
    7   3.8e+00  2.3e-02  1.9e-01  -5.79e-01  2.431154461e+04   2.432723725e+04   2.3e-02  104.31
    8   2.5e+00  1.5e-02  1.5e-01  -6.10e-01  2.732513220e+04   2.735001169e+04   1.5e-02  104.85
    9   8.2e-01  5.0e-03  4.6e-02  -2.83e-01  9.010287136e+04   9.012335414e+04   5.0e-03  105.39
    10  3.5e-01  2.1e-03  1.4e-02  3.53e-01   1.125156818e+05   1.125269640e+05   2.1e-03  105.94
    11  2.5e-02  1.5e-04  3.9e-04  6.50e-01   1.386274060e+05   1.386289381e+05   1.5e-04  106.70
    12  1.8e-02  1.1e-04  2.4e-04  9.63e-01   1.395665562e+05   1.395676758e+05   1.1e-04  107.24
    13  1.1e-03  7.0e-06  3.7e-06  9.76e-01   1.419509636e+05   1.419510327e+05   7.0e-06  107.83
    14  8.8e-06  5.4e-08  5.7e-09  9.98e-01   1.421098259e+05   1.421098265e+05   5.4e-08  108.90
    15  7.7e-06  4.7e-08  4.4e-09  1.00e+00   1.421099807e+05   1.421099812e+05   4.7e-08  109.88
    16  7.7e-06  4.7e-08  1.5e-09  1.00e+00   1.421099809e+05   1.421099814e+05   4.7e-08  111.03
    17  7.6e-06  3.5e-08  6.0e-10  9.99e-01   1.421102513e+05   1.421102516e+05   3.5e-08  112.09
    18  7.6e-06  3.5e-08  6.0e-10  1.00e+00   1.421102513e+05   1.421102516e+05   3.5e-08  113.34
    19  7.6e-06  3.5e-08  6.0e-10  1.00e+00   1.421102513e+05   1.421102516e+05   3.5e-08  114.66
    20  6.6e-06  3.1e-08  3.0e-10  1.00e+00   1.421103529e+05   1.421103532e+05   3.1e-08  115.64
    21  6.6e-06  3.1e-08  3.0e-10  1.00e+00   1.421103529e+05   1.421103532e+05   3.1e-08  116.84
    22  6.6e-06  3.1e-08  3.0e-10  1.00e+00   1.421103529e+05   1.421103532e+05   3.1e-08  118.15
    23  6.8e-06  3.0e-08  6.9e-11  1.00e+00   1.421103640e+05   1.421103643e+05   3.0e-08  119.25
    24  6.8e-06  3.0e-08  6.9e-11  1.00e+00   1.421103640e+05   1.421103643e+05   3.0e-08  120.44
    Optimizer terminated. Time: 122.08  


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.4211036402e+05    nrm: 3e+05    Viol.  con: 3e-06    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.4211036433e+05    nrm: 6e+05    Viol.  con: 2e-02    var: 1e-08    cones: 0e+00  



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

   **Total running time of the script:** ( 4 minutes  26.250 seconds)

**Estimated memory usage:**  4532 MB


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
