.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_coil_design_spherical_harmonics_coil_design.py>`     to download the full example code
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

    from bfieldtools.conductor import Conductor
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
    Computing the laplacian matrix...
    Computing the mass matrix...

    <mayavi.modules.vectors.Vectors object at 0x000001310E5183B8>



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
    Computing inductance matrix in 80 chunks (8840 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 56.65 seconds.
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
    Optimizer  - Scalar variables       : 110               conic                  : 82              
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.00              dense det. time        : 0.00            
    Factor     - ML order time          : 0.00              GP order time          : 0.00            
    Factor     - nonzeros before factor : 3321              after factor           : 3321            
    Factor     - dense dim.             : 0                 flops                  : 7.96e+05        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   4.0e+00  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  0.02  
    1   9.8e-01  2.4e-01  4.4e-01  -4.16e-02  1.125986870e+00   1.036083846e+00   2.4e-01  0.02  
    2   2.0e-01  5.0e-02  4.4e-02  3.29e-01   4.793669109e+00   4.753626457e+00   5.0e-02  0.02  
    3   9.4e-03  2.3e-03  1.2e-04  1.42e+00   5.783225382e+00   5.775510486e+00   2.3e-03  0.02  
    4   8.0e-04  2.0e-04  8.7e-06  1.22e+00   5.851111015e+00   5.850951896e+00   2.0e-04  0.02  
    5   1.2e-05  3.0e-06  1.6e-08  1.01e+00   5.852809333e+00   5.852806982e+00   3.0e-06  0.02  
    6   7.6e-08  1.9e-08  7.9e-12  1.00e+00   5.852851889e+00   5.852851874e+00   1.9e-08  0.02  
    7   4.2e-10  6.5e-11  2.0e-15  1.00e+00   5.852852242e+00   5.852852241e+00   8.7e-12  0.02  
    Optimizer terminated. Time: 0.02    


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 5.8528522417e+00    nrm: 1e+01    Viol.  con: 1e-11    var: 0e+00    cones: 0e+00  
      Dual.    obj: 5.8528522410e+00    nrm: 2e+02    Viol.  con: 2e-09    var: 1e-10    cones: 0e+00  




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

    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.33 seconds.

    <mayavi.modules.vectors.Vectors object at 0x000001310F7BB518>




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 1 minutes  23.821 seconds)


.. _sphx_glr_download_auto_examples_coil_design_spherical_harmonics_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: spherical_harmonics_coil_design.py <spherical_harmonics_coil_design.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: spherical_harmonics_coil_design.ipynb <spherical_harmonics_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
