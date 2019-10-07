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


    import numpy as np
    import matplotlib.pyplot as plt
    from mayavi import mlab
    import trimesh

    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.magnetic_field_mesh import compute_C
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops


    from bfieldtools.sphtools import compute_sphcoeffs_mesh, sphbasis, plotsph, sphfittools


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

    lmax = 10
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

    target_alms[0] += 1


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




    sph = sphbasis(10)
    sphfield = sph.field(target_points,target_alms, target_blms, lmax)

    target_field = sphfield/np.max(sphfield[:, 0])

    target_field[:, 2] = 0

    coil.plot_mesh()
    mlab.quiver3d(*target_points.T, *sphfield.T)






.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    l = 1 computed
    l = 2 computed
    l = 3 computed
    l = 4 computed
    l = 5 computed
    l = 6 computed
    l = 7 computed
    l = 8 computed
    l = 9 computed
    l = 10 computed



Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function

    #target_abs_error = np.zeros_like(target_alms)
    #target_abs_error += 0.5

    target_field = np.zeros_like(target_points)
    target_field[:, 0] += 1

    target_abs_error = np.zeros_like(target_points)
    target_abs_error += 0.01

    coil.C = compute_C(coil.mesh, target_points)

    #target_spec = {'C':coil.C_alms_norm, 'rel_error':None, 'abs_error':target_abs_error, 'target_field':target_alms}
    target_spec = {'C':coil.C, 'rel_error':None, 'abs_error':target_abs_error, 'target_field':target_field}






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 3184 vertices by 160 target points... took 0.30 seconds.



Run QP solver


.. code-block:: default

    import mosek

    coil.I, prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )



    B_target = coil.C.transpose([0, 2, 1]) @ coil.I


    lmax = 4
    coil.C_alms, coil.C_blms = compute_sphcoeffs_mesh(coil.mesh, lmax=lmax)

    Alms, Blms = coil.C_alms @ coil.I, coil.C_blms @ coil.I

    Blms = np.zeros_like(Alms)
    sphfield_target = sph.field(target_points, Alms, Blms, lmax)


    coeffs, coeffs2, nrmse = sphfittools.fitSpectra(sph, np.repeat(target_points[:, :, None], 3, -1), B_target, lmax)







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 2 chunks since 9 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 69.94 seconds.


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 3858            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 3858            
      Cones                  : 1               
      Scalar variables       : 5795            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 2897
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 3858              conic                  : 2898            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 1.06              dense det. time        : 0.00            
    Factor     - ML order time          : 0.28              GP order time          : 0.00            
    Factor     - nonzeros before factor : 4.20e+06          after factor           : 4.20e+06        
    Factor     - dense dim.             : 0                 flops                  : 3.64e+10        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   6.5e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  87.28 
    1   3.2e+01  5.0e-01  4.4e-01  2.67e-01   2.612520399e+01   2.542986715e+01   5.0e-01  87.69 
    2   1.0e+01  1.6e-01  5.0e-02  5.98e-01   6.044181048e+01   6.023015406e+01   1.6e-01  88.09 
    3   1.4e+00  2.1e-02  4.3e-03  1.15e+00   6.031796745e+01   6.029734437e+01   2.1e-02  88.57 
    4   5.9e-01  9.1e-03  1.2e-03  1.04e+00   6.115397807e+01   6.114526952e+01   9.1e-03  89.00 
    5   1.4e-02  2.1e-04  4.0e-06  1.01e+00   6.177515743e+01   6.177494195e+01   2.1e-04  89.46 
    6   6.4e-03  9.9e-05  1.3e-06  1.00e+00   6.178488084e+01   6.178478271e+01   9.9e-05  89.85 
    7   3.0e-03  4.7e-05  4.3e-07  1.00e+00   6.178656779e+01   6.178652232e+01   4.7e-05  90.23 
    8   3.4e-06  5.3e-08  1.3e-11  1.00e+00   6.178995680e+01   6.178995674e+01   5.3e-08  90.74 
    9   7.6e-07  2.4e-09  6.0e-13  1.00e+00   6.178996045e+01   6.178996031e+01   2.8e-11  91.37 
    Optimizer terminated. Time: 91.64   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 6.1789960450e+01    nrm: 1e+02    Viol.  con: 2e-11    var: 0e+00    cones: 0e+00  
      Dual.    obj: 6.1789960312e+01    nrm: 1e+02    Viol.  con: 1e-09    var: 4e-10    cones: 0e+00  
    l = 1 computed
    l = 2 computed
    l = 3 computed
    l = 4 computed
    Condition number = 1.890216
    Normalized RMS error = 0.216857%



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


.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_coil_design_002.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 4 minutes  12.881 seconds)

**Estimated memory usage:**  7831 MB


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
