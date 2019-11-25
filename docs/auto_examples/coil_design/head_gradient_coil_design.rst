.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_head_gradient_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_head_gradient_coil_design.py:


Head gradient coil
==================

Example showing a gradient coil designed on the surface of a MEG system helmet


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
    helmetmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools',
                                                                       'example_meshes/meg_helmet.obj'),
                              process=False)

    #planemesh.apply_scale(scaling_factor)
    #
    ##Specify coil plane geometry
    #center_offset = np.array([0, 0, 0]) * scaling_factor
    #standoff = np.array([0, 4, 0]) * scaling_factor
    #
    ##Create coil plane pairs
    #coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
    #                         planemesh.faces, process=False)
    ##
    #coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
    #                     planemesh.faces, process=False)

    #joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = MeshWrapper(verts=helmetmesh.vertices, tris=helmetmesh.faces, fix_normals=True)







Set up target and stray field points.
Here, the target points are on a volumetric grid within a sphere


.. code-block:: default


    offset = np.array([0, 0, 0.04])
    center = offset * scaling_factor

    sidelength = 0.05 * scaling_factor
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









Specify target field and run solver


.. code-block:: default


    #Let's generate the target field through the use of spherical harmonics.
    # Thus we avoid issues with having to manually specify the concomitant gradients


    from bfieldtools.sphtools import sphbasis


    sph = sphbasis(50)

    #plotsph.plotYlms(sph, 3)

    lmax = 3
    alm = np.zeros((lmax*(lmax+2),))
    blm = np.zeros((lmax*(lmax+2),))

    #

    blm[3]+=1

    sphfield = sph.field(target_points - offset, alm, blm, lmax)

    target_field = sphfield/np.max(sphfield[:, 0])

    target_field[:, 2] = 0

    coil.plot_mesh()
    mlab.quiver3d(*target_points.T, *sphfield.T)



    rel_error = np.zeros_like(target_field)
    #rel_error[:, 0] += 0.1

    abs_error = np.zeros_like(target_field)
    abs_error[:, 0] += 0.1
    abs_error[:, 1:3] += 0.1


    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':rel_error, 'abs_error':abs_error, 'target':target_field}

    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )




.. image:: /auto_examples/coil_design/images/sphx_glr_head_gradient_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 2044 vertices by 672 target points... took 0.62 seconds.
    Computing inductance matrix in 2 chunks since 6 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 31.69 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 5979            
      Cones                  : 1               
      Scalar variables       : 3893            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 5979            
      Cones                  : 1               
      Scalar variables       : 3893            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 1946
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 5979              conic                  : 1947            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.97              dense det. time        : 0.00            
    Factor     - ML order time          : 0.10              GP order time          : 0.00            
    Factor     - nonzeros before factor : 1.89e+06          after factor           : 1.89e+06        
    Factor     - dense dim.             : 0                 flops                  : 1.75e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   2.8e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  39.06 
    1   9.6e+01  3.4e-01  1.0e+00  -8.89e-01  3.495137057e+01   3.556319782e+01   3.4e-01  39.39 
    2   1.9e+01  6.7e-02  3.1e-01  -6.80e-01  3.287714246e+02   3.330588000e+02   6.7e-02  39.66 
    3   3.0e+00  1.1e-02  3.5e-02  -3.87e-03  7.839861012e+02   7.864354333e+02   1.1e-02  39.94 
    4   2.4e-01  8.5e-04  7.5e-04  7.33e-01   7.252185038e+02   7.253864548e+02   8.5e-04  40.31 
    5   2.2e-01  7.9e-04  6.8e-04  9.90e-01   7.113676359e+02   7.115249096e+02   7.9e-04  40.57 
    6   2.1e-01  7.6e-04  6.3e-04  9.87e-01   7.125757105e+02   7.127264392e+02   7.6e-04  40.84 
    7   1.4e-01  4.9e-04  3.4e-04  9.86e-01   6.702376812e+02   6.703372912e+02   4.9e-04  41.07 
    8   9.2e-02  3.3e-04  1.8e-04  9.89e-01   6.598383453e+02   6.599050758e+02   3.3e-04  41.37 
    9   5.4e-02  1.9e-04  8.2e-05  9.92e-01   6.563197990e+02   6.563594702e+02   1.9e-04  41.63 
    10  2.5e-02  8.8e-05  2.5e-05  9.95e-01   6.509659207e+02   6.509842545e+02   8.8e-05  41.90 
    11  1.4e-02  5.1e-05  1.1e-05  9.98e-01   6.498164921e+02   6.498271837e+02   5.1e-05  42.14 
    12  1.6e-03  5.8e-06  4.4e-07  9.99e-01   6.485899413e+02   6.485912086e+02   5.8e-06  42.46 
    13  3.4e-05  1.2e-07  1.3e-09  1.00e+00   6.486172624e+02   6.486172882e+02   1.2e-07  42.81 
    14  6.0e-06  2.1e-08  9.0e-11  1.00e+00   6.486175937e+02   6.486175984e+02   2.1e-08  43.05 
    15  2.5e-07  6.9e-10  1.4e-11  1.00e+00   6.486178409e+02   6.486178382e+02   2.6e-11  43.56 
    Optimizer terminated. Time: 43.82   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 6.4861784095e+02    nrm: 1e+03    Viol.  con: 5e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 6.4861783820e+02    nrm: 3e+03    Viol.  con: 7e-08    var: 5e-09    cones: 0e+00  



Plot coil windings and magnetic field in target points


.. code-block:: default



    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=20)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.05/50)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)

    f.scene.isometric_view()



.. image:: /auto_examples/coil_design/images/sphx_glr_head_gradient_coil_design_002.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 1 minutes  37.063 seconds)

**Estimated memory usage:**  3479 MB


.. _sphx_glr_download_auto_examples_coil_design_head_gradient_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: head_gradient_coil_design.py <head_gradient_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: head_gradient_coil_design.ipynb <head_gradient_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
