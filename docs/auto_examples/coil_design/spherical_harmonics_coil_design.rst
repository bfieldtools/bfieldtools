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

    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops


    from bfieldtools.sphtools import compute_sphcoeffs_mesh, sphbasis, sphfittools


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

    lmax = 4
    coil.C_alms, coil.C_blms = compute_sphcoeffs_mesh(coil.mesh, lmax=lmax)


    #Radius of sphere of interest
    Rmax = 1.0

    lind = 0
    coil.C_alms_norm = np.zeros_like(coil.C_alms)
    for l in range(1,lmax+1):
        for m in range(-1*l,l+1):
            temp = (2*l**2 + l)*Rmax**(2*l-1)/(2*l-1)
            #coeffs2[lind] = coeffs[lind]**2*temp
            coil.C_alms_norm[lind] = coil.C_alms[lind]/temp**0.5
            lind += 1

    target_alms = np.zeros((lmax * (lmax+2),))
    target_blms = np.zeros((lmax * (lmax+2),))

    target_blms[0] += 1


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




    sph = sphbasis(4)
    sphfield = sph.field(target_points, target_alms, target_blms, lmax)

    target_field = sphfield/np.max(sphfield[:, 0])

    target_field[:, 2] = 0

    coil.plot_mesh()
    mlab.quiver3d(*target_points.T, *sphfield.T)






.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    SVG path loading unavailable!
    Traceback (most recent call last):
      File "/u/76/zetterr1/unix/.local/lib/python3.6/site-packages/trimesh/path/exchange/svg_io.py", line 18, in <module>
        from svg.path import parse_path
    ModuleNotFoundError: No module named 'svg'
    l = 1 computed
    l = 2 computed
    l = 3 computed
    l = 4 computed



Create bfield specifications used when optimizing the coil geometry


.. code-block:: default



    target_spec = {'coupling':coil.C_blms, 'rel_error':0, 'abs_error':0.01, 'target':target_blms}








Run QP solver


.. code-block:: default

    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )

    B_target = coil.B_coupling(target_points) @ coil.j


    lmax = 4
    coil.C_alms, coil.C_blms = compute_sphcoeffs_mesh(coil.mesh, lmax=lmax)

    Alms, Blms = coil.C_alms @ coil.j, coil.C_blms @ coil.j

    Alms = np.zeros_like(Blms)
    sphfield_target = sph.field(target_points, Alms, Blms, lmax)


    coeffs, coeffs2, nrmse = sphfittools.fitSpectra(sph, np.repeat(target_points[:, :, None], 3, -1), B_target, lmax)







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 2 chunks since 8 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 70.08 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 2946            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 2946            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 2946              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.90              dense det. time        : 0.00            
    Factor     - ML order time          : 0.30              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 3.26e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.0e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  40.05 
    1   1.4e+00  1.4e-01  2.9e-01  -5.70e-01  4.541786587e-01   9.355049528e-01   1.4e-01  40.52 
    2   2.6e-01  2.6e-02  1.5e-02  1.01e+00   3.792613273e-01   3.617956708e-01   2.6e-02  40.87 
    3   1.1e-02  1.1e-03  9.4e-05  1.20e+00   3.063693271e-01   3.040294637e-01   1.1e-03  41.32 
    4   2.5e-04  2.4e-05  3.8e-07  1.04e+00   2.902837007e-01   2.902580349e-01   2.4e-05  41.76 
    5   9.4e-06  9.1e-07  3.1e-09  1.00e+00   2.894031176e-01   2.894027587e-01   9.1e-07  42.14 
    6   5.5e-08  5.4e-09  1.4e-12  1.00e+00   2.893970282e-01   2.893970253e-01   5.4e-09  42.64 
    Optimizer terminated. Time: 42.85   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 2.8939702817e-01    nrm: 2e+00    Viol.  con: 1e-08    var: 0e+00    cones: 0e+00  
      Dual.    obj: 2.8939702529e-01    nrm: 1e+01    Viol.  con: 8e-06    var: 6e-14    cones: 0e+00  
    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.25 seconds.
    l = 1 computed
    l = 2 computed
    l = 3 computed
    l = 4 computed
    Condition number = 1.890216
    Normalized RMS error = 0.157768%



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


.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_coil_design_002.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  57.785 seconds)

**Estimated memory usage:**  7964 MB


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
