.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_self-shielded_biplanar_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_self-shielded_biplanar_coil_design.py:


Analytical self-shielded biplanar coil design
==============================================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes. In addition, the coils have an outer surface
for which (in a linear fashion) a secondary current is created, which zeroes the
normal component of the field produced by the primary coil at the secondary coil
surface. The combination of the primary and secondary coil currents are specified to create
the target field, and their combined inductive energy is minimized.

NB. The secondary coil current is entirely a function of the primary coil current
and the geometry.


.. code-block:: default


    import numpy as np
    import matplotlib.pyplot as plt
    from mayavi import mlab
    import trimesh

    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic, scalar_potential_coupling
    from bfieldtools.mesh_properties import mutual_inductance_matrix
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops, plot_data_on_vertices

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane.obj'), process=False)

    planemesh.apply_scale(scaling_factor*1.6)

    #Specify coil plane geometry
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 5, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                             planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                         planemesh.faces, process=False)

    joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = MeshWrapper(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)

    shieldmesh = joined_planes.copy()
    shieldmesh.vertices *= np.array([1.5, 1.5, 1.5])

    shieldcoil = MeshWrapper(verts=shieldmesh.vertices, tris=shieldmesh.faces, fix_normals=True)







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    SVG path loading unavailable!
    Traceback (most recent call last):
      File "/u/76/zetterr1/unix/.local/lib/python3.6/site-packages/trimesh/path/exchange/svg_io.py", line 18, in <module>
        from svg.path import parse_path
    ModuleNotFoundError: No module named 'svg'



Compute inductances and coupling


.. code-block:: default



    M11 = coil.inductance
    M22 = shieldcoil.inductance
    # Constrain boundary to zero and consider only inneverts
    M11 = M11#[coil.inner_verts][:, coil.inner_verts]
    M22 = M22[shieldcoil.inner_verts][:, shieldcoil.inner_verts]
    # Add rank-one matrix, so that M22 can be inverted (for zero mean functions)
    #M22 += np.ones_like(M22)/M22.shape[0]
    #M11 += np.ones_like(M11)/M11.shape[0]



    M21 = mutual_inductance_matrix(shieldcoil.mesh, coil.mesh)
    M21 = M21[shieldcoil.inner_verts]

    # Mapping from I1 to I2, constraining flux through shieldcoil to zero
    P = -np.linalg.solve(M22, M21)







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 73116 MiB required for 1352 times 1352 vertices...
    Computing inductance matrix in 8 chunks since 9993 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 14.42 seconds.
    Computing self-inductance matrix using rough quadrature. For higher accuracy, set quad_degree to 4 or more.
    Estimating 73116 MiB required for 1352 times 1352 vertices...
    Computing inductance matrix in 8 chunks since 9944 MiB memory is available...
    Computing potential matrix
    Inductance matrix computation took 14.21 seconds.
    Estimating 73116 MiB required for 1352 times 1352 vertices...
    Computing inductance matrix in 8 chunks since 9915 MiB memory is available...
    Computing potential matrix



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








Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function

    target_field = np.zeros(target_points.shape)
    target_field[:, 1] = target_field[:, 1] + 1

    target_rel_error = np.zeros_like(target_field)
    target_rel_error[:, 0] += 0.01

    target_abs_error = np.zeros_like(target_field)
    target_abs_error[:, 0] += 0.001
    target_abs_error[:, 1:3] += 0.005

    target_spec = {'coupling':coil.B_coupling(target_points) + shieldcoil.B_coupling(target_points)[:, :, shieldcoil.inner_verts]@P, 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}
    #[:, :, coil.inner_verts]

    objective_matrix = M11 - M21.T @ np.linalg.pinv(M22) @ M21





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 1352 vertices by 160 target points... took 0.12 seconds.
    Computing magnetic field coupling matrix, 1352 vertices by 160 target points... took 0.11 seconds.



Run QP solver


.. code-block:: default

    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective=objective_matrix,
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}},
                                       boundary_constraints='all_zero'
                                       )

    shieldcoil.j = np.zeros((len(shieldcoil.mesh.vertices, )))

    shieldcoil.j[shieldcoil.inner_verts] = P @ coil.j



    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))

    plot_data_on_vertices(coil.mesh, coil.j, figure=f)
    plot_data_on_vertices(shieldcoil.mesh, shieldcoil.j, figure=f)




.. image:: /auto_examples/coil_design/images/sphx_glr_self-shielded_biplanar_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /l/bfieldtools/bfieldtools/coil_optimize.py:175: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
      if objective == 'minimum_inductive_energy':
    /l/bfieldtools/bfieldtools/coil_optimize.py:177: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
      elif objective == 'minimum_resistive_energy':
    Custom objective passed, assuming it is a matrix of correct dimensions
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...
    /l/conda-envs/mne/lib/python3.6/site-packages/cvxpy/reductions/solvers/solving_chain.py:170: UserWarning: You are solving a parameterized problem that is not DPP. Because the problem is not DPP, subsequent solves will not be faster than the first one.
      "You are solving a parameterized problem that is not DPP. "


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 2130            
      Cones                  : 1               
      Scalar variables       : 2339            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 2130            
      Cones                  : 1               
      Scalar variables       : 2339            
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 1169
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 2130              conic                  : 1170            
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.13              dense det. time        : 0.00            
    Factor     - ML order time          : 0.02              GP order time          : 0.00            
    Factor     - nonzeros before factor : 6.84e+05          after factor           : 6.84e+05        
    Factor     - dense dim.             : 0                 flops                  : 2.78e+09        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   1.6e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  7.89  
    1   1.2e+01  7.3e-01  1.0e+00  1.87e-01   3.204310871e+00   2.373758489e+00   7.3e-01  7.95  
    2   8.6e+00  5.4e-01  1.3e-01  4.10e-01   1.173977597e+01   1.117392378e+01   5.4e-01  8.01  
    3   6.8e+00  4.2e-01  1.1e-01  3.65e+00   1.683957085e+01   1.653068413e+01   4.2e-01  8.06  
    4   2.8e+00  1.7e-01  3.9e-02  2.86e+00   1.804535013e+01   1.797918789e+01   1.7e-01  8.11  
    5   1.6e+00  9.7e-02  1.5e-02  1.53e+00   1.849297968e+01   1.846055406e+01   9.7e-02  8.16  
    6   8.8e-02  5.4e-03  1.7e-04  1.32e+00   1.885791271e+01   1.885630498e+01   5.4e-03  8.24  
    7   5.6e-03  3.5e-04  2.2e-06  1.03e+00   1.889865051e+01   1.889854392e+01   3.5e-04  8.30  
    8   3.1e-06  1.9e-07  2.9e-11  1.00e+00   1.890146314e+01   1.890146308e+01   1.9e-07  8.38  
    9   9.8e-08  1.5e-09  1.9e-13  1.00e+00   1.890146467e+01   1.890146466e+01   6.6e-10  8.50  
    Optimizer terminated. Time: 8.55    


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 1.8901464672e+01    nrm: 4e+01    Viol.  con: 1e-10    var: 0e+00    cones: 0e+00  
      Dual.    obj: 1.8901464665e+01    nrm: 5e+01    Viol.  con: 2e-09    var: 6e-11    cones: 0e+00  



Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=N_contours)
    sloops, sloop_values= scalar_contour(shieldcoil.mesh, shieldcoil.j, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)
    plot_3d_current_loops(sloops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.j + shieldcoil.B_coupling(target_points) @ shieldcoil.j

    mlab.quiver3d(*target_points.T, *B_target.T)




    extent = 30




.. image:: /auto_examples/coil_design/images/sphx_glr_self-shielded_biplanar_coil_design_002.png
    :class: sphx-glr-single-img




Compute field along major axes


.. code-block:: default



    x1 = np.linspace(-extent, extent, 101) * scaling_factor

    y1 = z1 = np.zeros_like(x1)

    line1_points = np.vstack((x1, y1, z1)).T

    B_line1 = coil.B_coupling(line1_points) @ coil.j + shieldcoil.B_coupling(line1_points) @ shieldcoil.j


    y2 = np.linspace(-extent, extent, 101) * scaling_factor

    z2 = x2 = np.zeros_like(y2)

    line2_points = np.vstack((x2, y2, z2)).T

    B_line2 = coil.B_coupling(line2_points) @ coil.j + shieldcoil.B_coupling(line2_points) @ shieldcoil.j



    z3 = np.linspace(-extent, extent, 101) * scaling_factor

    x3 = y3 = np.zeros_like(z1)

    line3_points = np.vstack((x3, y3, z3)).T


    B_line3 = coil.B_coupling(line3_points) @ coil.j + shieldcoil.B_coupling(line3_points) @ shieldcoil.j

    fig, axes = plt.subplots(1, 1)

    for ax_idx, ax in enumerate([axes]):
        ax.semilogy(x1 / scaling_factor, np.linalg.norm(B_line1, axis=-1), label='X')
        ax.semilogy(y2 / scaling_factor, np.linalg.norm(B_line2, axis=-1), label='Y')
        ax.semilogy(z3 / scaling_factor, np.linalg.norm(B_line3, axis=-1), label='Z')
        ax.set_title('Field component %d'% ax_idx)

    plt.ylabel('Field amplitude (target field units)')
    plt.xlabel('Distance from origin')
    plt.grid(True, which='minor', axis='y')
    plt.grid(True, which='major', axis='y', color='k')
    plt.grid(True, which='major', axis='x')

    plt.legend()


    plt.show()




.. image:: /auto_examples/coil_design/images/sphx_glr_self-shielded_biplanar_coil_design_003.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 1352 vertices by 101 target points... took 0.09 seconds.
    Computing magnetic field coupling matrix, 1352 vertices by 101 target points... took 0.08 seconds.
    Computing magnetic field coupling matrix, 1352 vertices by 100 target points... took 0.08 seconds.
    Computing magnetic field coupling matrix, 1352 vertices by 100 target points... took 0.08 seconds.
    Computing magnetic field coupling matrix, 1352 vertices by 100 target points... took 0.08 seconds.
    Computing magnetic field coupling matrix, 1352 vertices by 100 target points... took 0.08 seconds.
    /l/bfieldtools/examples/coil_design/self-shielded_biplanar_coil_design.py:229: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
      plt.show()



Compute the field and scalar potential on a larger plane


.. code-block:: default


    x = y = np.linspace(-20, 20, 50)
    X,Y = np.meshgrid(x, y, indexing='ij')
    points = np.zeros((X.flatten().shape[0], 3))
    points[:, 0] = X.flatten()
    points[:, 1] = Y.flatten()

    CB1 = magnetic_field_coupling_analytic(coil.mesh, points)
    CB2 = magnetic_field_coupling_analytic(shieldcoil.mesh, points)

    CU1 = scalar_potential_coupling(coil.mesh, points)
    CU2 = scalar_potential_coupling(shieldcoil.mesh, points)

    B1 = CB1 @ coil.j
    B2 = CB2 @ shieldcoil.j

    U1 = CU1 @ coil.j
    U2 = CU2 @ shieldcoil.j






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix analytically, 1352 vertices by 2500 target points... took 4.73 seconds.
    Computing magnetic field coupling matrix analytically, 1352 vertices by 2500 target points... took 4.76 seconds.
    Computing scalar potential coupling matrix, 1352 vertices by 2500 target points... took 6.32 seconds.
    Computing scalar potential coupling matrix, 1352 vertices by 2500 target points... took 6.26 seconds.



Plot field and potential planar cross-section


.. code-block:: default

    B = (B1.T + B2.T)[:2].reshape(2, x.shape[0], y.shape[0])
    lw = np.sqrt(B[0]**2 + B[1]**2)
    lw = 2*lw/np.max(lw)
    xx = np.linspace(-1,1, 16)
    #seed_points = 0.51*np.array([xx, -np.sqrt(1-xx**2)])
    #seed_points = np.hstack([seed_points, (0.51*np.array([xx, np.sqrt(1-xx**2)]))])
    #plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
    #               start_points=seed_points.T, integration_direction='both')
    U = (U1 + U2).reshape(x.shape[0], y.shape[0])
    U /= np.max(U)
    plt.figure()
    plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
               extent=(x.min(), x.max(), y.min(), y.max()))
    plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
                   #start_points=seed_points.T,
                   integration_direction='both')

    cc1 = scalar_contour(coil.mesh, coil.mesh.vertices[:,2], contours= [-0.001])[0][0]
    cc2 = scalar_contour(shieldcoil.mesh, shieldcoil.mesh.vertices[:,2], contours= [-0.001])[0][0]

    plt.plot(cc1[:,1], cc1[:,0], linewidth=3.0)
    plt.plot(cc2[:,1], cc2[:,0], linewidth=3.0)

    plt.xticks([])
    plt.yticks([])






.. image:: /auto_examples/coil_design/images/sphx_glr_self-shielded_biplanar_coil_design_004.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 1 minutes  16.750 seconds)

**Estimated memory usage:**  1113 MB


.. _sphx_glr_download_auto_examples_coil_design_self-shielded_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: self-shielded_biplanar_coil_design.py <self-shielded_biplanar_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: self-shielded_biplanar_coil_design.ipynb <self-shielded_biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
