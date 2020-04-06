.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_coil_design_head_gradient_coil_design.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_coil_design_head_gradient_coil_design.py:


Head gradient coil
==================

Example showing a gradient coil designed on the surface of a MEG system helmet


.. code-block:: default



    import numpy as np
    from mayavi import mlab
    import trimesh


    from bfieldtools.conductor import Conductor
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
    coil = Conductor(verts=helmetmesh.vertices, tris=helmetmesh.faces, fix_normals=True)








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


    from bfieldtools import sphtools


    lmax = 3
    alm = np.zeros((lmax*(lmax+2),))
    blm = np.zeros((lmax*(lmax+2),))

    #

    blm[3]+=1

    sphfield = sphtools.field(target_points - offset, alm, blm, lmax)

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

    Computing magnetic field coupling matrix, 2044 vertices by 672 target points... took 0.75 seconds.
    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 16313 MiB required for 2044 by 2044 vertices...
    Computing inductance matrix in 80 chunks (5418 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 25.88 seconds.
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
    Factor     - setup time             : 0.58              dense det. time        : 0.00            
    Factor     - ML order time          : 0.09              GP order time          : 0.00            
    Factor     - nonzeros before factor : 1.89e+06          after factor           : 1.89e+06        
    Factor     - dense dim.             : 0                 flops                  : 1.75e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.4e+02  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  50.16 
    1   5.1e+01  3.6e-01  1.0e+00  -8.22e-01  5.753343851e+01   5.786493454e+01   3.6e-01  50.72 
    2   1.4e+01  9.7e-02  3.4e-01  -5.55e-01  4.500836472e+02   4.522180008e+02   9.7e-02  51.22 
    3   1.8e+00  1.3e-02  3.2e-02  5.33e-02   9.768646705e+02   9.781318367e+02   1.3e-02  51.73 
    4   3.1e-01  2.2e-03  2.0e-03  8.17e-01   9.221917278e+02   9.223518215e+02   2.2e-03  52.34 
    5   8.9e-02  6.3e-04  3.1e-04  9.69e-01   9.162883526e+02   9.163346487e+02   6.3e-04  52.91 
    6   4.1e-02  2.9e-04  9.7e-05  9.90e-01   9.159409663e+02   9.159627834e+02   2.9e-04  53.45 
    7   2.4e-02  1.7e-04  4.4e-05  9.95e-01   9.165772787e+02   9.165901906e+02   1.7e-04  53.95 
    8   1.4e-03  9.6e-06  5.9e-07  9.97e-01   9.173529626e+02   9.173537201e+02   9.6e-06  54.53 
    9   2.3e-05  1.7e-07  1.3e-09  1.00e+00   9.174520496e+02   9.174520627e+02   1.7e-07  55.08 
    10  8.0e-07  1.9e-09  7.9e-12  1.00e+00   9.174541177e+02   9.174541151e+02   9.6e-12  55.61 
    Optimizer terminated. Time: 55.91   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 9.1745411770e+02    nrm: 2e+03    Viol.  con: 1e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 9.1745411514e+02    nrm: 4e+03    Viol.  con: 2e-08    var: 4e-09    cones: 0e+00  




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

   **Total running time of the script:** ( 1 minutes  50.616 seconds)


.. _sphx_glr_download_auto_examples_coil_design_head_gradient_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: head_gradient_coil_design.py <head_gradient_coil_design.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: head_gradient_coil_design.ipynb <head_gradient_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
