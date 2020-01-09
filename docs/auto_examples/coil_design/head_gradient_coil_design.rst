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

    Computing magnetic field coupling matrix, 2044 vertices by 672 target points... took 0.48 seconds.
    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 167117 MiB required for 2044 times 2044 vertices...
    Computing inductance matrix in 17 chunks since 10116 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 35.78 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


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
    Factor     - setup time             : 0.94              dense det. time        : 0.00            
    Factor     - ML order time          : 0.10              GP order time          : 0.00            
    Factor     - nonzeros before factor : 1.89e+06          after factor           : 1.89e+06        
    Factor     - dense dim.             : 0                 flops                  : 1.75e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   2.8e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  43.22 
    1   1.1e+02  3.8e-01  1.1e+00  -8.59e-01  4.052113853e+01   4.085107520e+01   3.8e-01  43.48 
    2   1.6e+01  5.8e-02  2.6e-01  -6.52e-01  4.812095471e+02   4.851725448e+02   5.8e-02  43.73 
    3   2.8e+00  9.8e-03  2.8e-02  1.76e-01   8.764755952e+02   8.782528059e+02   9.8e-03  43.96 
    4   7.2e-01  2.6e-03  3.8e-03  8.05e-01   7.942959523e+02   7.947842332e+02   2.6e-03  44.18 
    5   4.8e-01  1.7e-03  2.1e-03  9.50e-01   7.682595053e+02   7.685747069e+02   1.7e-03  44.41 
    6   8.3e-02  3.0e-04  1.3e-04  9.65e-01   7.541048948e+02   7.541438973e+02   3.0e-04  44.72 
    7   4.3e-02  1.5e-04  4.7e-05  9.93e-01   7.541445121e+02   7.541646991e+02   1.5e-04  44.94 
    8   7.2e-03  2.6e-05  3.3e-06  9.97e-01   7.546455334e+02   7.546491101e+02   2.6e-05  45.19 
    9   9.5e-04  3.4e-06  1.6e-07  9.99e-01   7.548605518e+02   7.548610249e+02   3.4e-06  45.46 
    10  3.3e-04  1.2e-06  3.2e-08  1.00e+00   7.548815988e+02   7.548817620e+02   1.2e-06  45.71 
    11  3.8e-05  1.3e-07  1.3e-09  1.00e+00   7.548926602e+02   7.548926790e+02   1.3e-07  45.95 
    12  3.7e-06  1.3e-08  4.2e-11  1.00e+00   7.548941296e+02   7.548941314e+02   1.3e-08  46.45 
    13  5.0e-07  1.8e-09  3.8e-11  1.00e+00   7.548942714e+02   7.548942725e+02   1.8e-09  46.67 
    14  2.7e-07  9.4e-10  5.7e-11  1.00e+00   7.548942818e+02   7.548942801e+02   9.4e-10  47.02 
    15  2.5e-07  4.9e-10  2.8e-11  1.00e+00   7.548942875e+02   7.548942863e+02   4.9e-10  47.41 
    Optimizer terminated. Time: 47.68   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 7.5489428748e+02    nrm: 2e+03    Viol.  con: 7e-09    var: 0e+00    cones: 0e+00  
      Dual.    obj: 7.5489428632e+02    nrm: 4e+03    Viol.  con: 1e-06    var: 2e-09    cones: 0e+00  



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

   **Total running time of the script:** ( 1 minutes  45.173 seconds)

**Estimated memory usage:**  2239 MB


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
