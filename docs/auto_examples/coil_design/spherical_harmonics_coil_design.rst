.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_spherical_harmonics_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_spherical_harmonics_coil_design.py:


Spherical harmonics-generating coil design
==========================================

Example showing a basic biplanar coil producing a field profile defined by
spherical harmonics.



.. code-block:: default


    #import sys
    #path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
    #if path in sys.path:
    #    sys.path.insert(0, path)


    import numpy as np
    from mayavi import mlab
    import trimesh

    from bfieldtools.mesh_class import Conductor
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops
    from bfieldtools import sphtools


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

    coil = Conductor(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True, N_sph=4, basis_name='suh', N_suh=80)

    target_alms = np.zeros((coil.opts['N_sph'] * (coil.opts['N_sph']+2),))
    target_blms = np.zeros((coil.opts['N_sph'] * (coil.opts['N_sph']+2),))

    target_blms[3] += 1


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



    sphfield = sphtools.field(target_points, target_alms, target_blms, coil.opts['N_sph'])


    target_field = sphfield/np.max(np.abs(sphfield))

    coil.plot_mesh()
    mlab.quiver3d(*target_points.T, *sphfield.T)






.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Calculating surface harmonics expansion...



Create bfield specifications used when optimizing the coil geometry


.. code-block:: default



    target_spec = {'coupling':coil.sph_couplings[1], 'rel_error':0, 'abs_error':0.01, 'target':target_blms}






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing coupling matrices
    l = 1 computed
    l = 2 computed
    l = 3 computed
    l = 4 computed



Run QP solver


.. code-block:: default

    import mosek

    coil.s, prob = optimize_streamfunctions(coil,
                                       [target_spec],
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
    Computing inductance matrix in 100 chunks (7503 MiB memory free), when approx_far=True using more chunks is faster...
    Computing potential matrix
    Inductance matrix computation took 52.38 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 130             
      Cones                  : 1               
      Scalar variables       : 163             
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 130             
      Cones                  : 1               
      Scalar variables       : 163             
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 81
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 130               conic                  : 82              
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.00              dense det. time        : 0.00            
    Factor     - ML order time          : 0.00              GP order time          : 0.00            
    Factor     - nonzeros before factor : 3321              after factor           : 3321            
    Factor     - dense dim.             : 0                 flops                  : 8.61e+05        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.0e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  0.03  
    1   3.6e+00  3.5e-01  7.8e-01  -6.16e-01  6.843060430e-01   5.766876460e-01   3.5e-01  0.07  
    2   1.7e+00  1.7e-01  3.8e-01  -2.47e-01  4.307239132e+00   4.600847171e+00   1.7e-01  0.07  
    3   7.7e-01  7.5e-02  1.4e-01  6.35e-02   1.122605585e+01   1.148745057e+01   7.5e-02  0.07  
    4   1.0e-01  1.0e-02  4.8e-03  9.38e-01   1.770154008e+01   1.768332845e+01   1.0e-02  0.07  
    5   1.5e-02  1.5e-03  2.2e-04  1.22e+00   1.883522060e+01   1.883138500e+01   1.5e-03  0.07  
    6   4.1e-04  4.0e-05  1.4e-06  1.12e+00   1.904763480e+01   1.904768955e+01   4.0e-05  0.07  
    7   6.9e-07  6.7e-08  9.3e-11  1.01e+00   1.905276377e+01   1.905276386e+01   6.7e-08  0.07  
    8   9.5e-10  9.5e-11  5.4e-15  1.00e+00   1.905277126e+01   1.905277127e+01   4.5e-12  0.08  
    Optimizer terminated. Time: 0.08    


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.9052771265e+01    nrm: 4e+01    Viol.  con: 1e-11    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.9052771268e+01    nrm: 6e+02    Viol.  con: 1e-08    var: 2e-10    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.s.vert, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.s

    mlab.quiver3d(*target_points.T, *B_target.T)


.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_coil_design_002.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.26 seconds.




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 1 minutes  12.342 seconds)


.. _sphx_glr_download_auto_examples_coil_design_spherical_harmonics_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: spherical_harmonics_coil_design.py <spherical_harmonics_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: spherical_harmonics_coil_design.ipynb <spherical_harmonics_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
