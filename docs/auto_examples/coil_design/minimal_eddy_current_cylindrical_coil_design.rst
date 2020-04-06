.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_coil_design_minimal_eddy_current_cylindrical_coil_design.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_coil_design_minimal_eddy_current_cylindrical_coil_design.py:


Coil with minimal eddy currents
===============================
Compact example of design of a cylindrical coil surrounded by a RF shield, i.e. a conductive surface.
The effects of eddy currents due to inductive interaction with the shield is minimized


.. code-block:: default


    import numpy as np
    from mayavi import mlab
    import trimesh


    from bfieldtools.conductor import Conductor

    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops, plot_data_on_vertices

    import pkg_resources

    from pyface.api import GUI
    _gui = GUI()


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load example coil mesh that is centered on the origin
    coilmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/open_cylinder.stl'), process=True)

    angle = np.pi/2
    rotation_matrix = np.array([[np.cos(angle), 0, np.sin(angle), 0],
                                  [0, 1, 0, 0],
                                  [-np.sin(angle), 0, np.cos(angle), 0],
                                  [0, 0, 0, 1]
                                  ])

    coilmesh.apply_transform(rotation_matrix)

    coilmesh1 = coilmesh.copy()
    #coilmesh1.apply_scale(1.3)

    coilmesh2 = coilmesh.copy()

    #coilmesh1 = coilmesh.union(coilmesh1)
    #coilmesh1 = coilmesh1.subdivide().subdivide()
    #coilmesh2 = coilmesh.subdivide()


    #Create mesh class object
    coil = Conductor(verts=coilmesh1.vertices*0.75, tris=coilmesh1.faces,
                     fix_normals=True,
                     basis_name='suh',
                     N_suh=400
                     )


    def alu_sigma(T):
        ref_T = 293 #K
        ref_rho = 2.82e-8 #ohm*meter
        alpha = 0.0039 #1/K


        rho = alpha * (T - ref_T) * ref_rho + ref_rho

        return 1/rho

    resistivity = 1/alu_sigma(T=293) #room-temp Aluminium
    thickness = 0.5e-3 # 0.5 mm thick


    # Separate object for shield geometry
    shield = Conductor(verts=coilmesh2.vertices.copy()*1.1, tris=coilmesh2.faces.copy(),
                       fix_normals=True,
                       basis_name='inner',
                       resistivity=resistivity,
                       thickness=thickness)







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    face_normals didn't match triangles, ignoring!
    Calculating surface harmonics expansion...
    Computing the laplacian matrix...
    Computing the mass matrix...




Set up target  points and plot geometry


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 0])

    sidelength = 0.25 * scaling_factor
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


    #Plot coil, shield and target points
    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                size=(800, 800))
    coil.plot_mesh()
    shield.plot_mesh()
    mlab.points3d(*target_points.T)







.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_001.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_003.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.glyph.Glyph object at 0x000001D585B1A150>



Compute C matrices that are used to compute the generated magnetic field


.. code-block:: default


    mutual_inductance = coil.mutual_inductance(shield)

    # Take into account the field produced by currents induced into the shield
    # NB! This expression is for instantaneous step-function switching of coil current, see Eq. 18 in G.N. Peeren, 2003.

    shield.M_coupling = np.linalg.solve(-shield.inductance, mutual_inductance.T)
    secondary_C = shield.B_coupling(target_points) @ -shield.M_coupling





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Estimating 69923 MiB required for 4764 by 4764 vertices...
    Computing inductance matrix in 180 chunks (8210 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 69923 MiB required for 4764 by 4764 vertices...
    Computing inductance matrix in 180 chunks (8305 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 126.19 seconds.
    Computing magnetic field coupling matrix, 4764 vertices by 672 target points... took 1.50 seconds.




Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function

    target_field = np.zeros(target_points.shape)
    target_field[:, 1] = target_field[:, 1] + 1


    target_spec = {'coupling':coil.B_coupling(target_points), 'abs_error':0.01, 'target':target_field}




    from scipy.linalg import eigh
    l, U = eigh(shield.resistance, shield.inductance, eigvals=(0, 500))
    #
    #U = np.zeros((shield.inductance.shape[0], len(li)))
    #U[shield.inner_verts, :] = Ui


    #
    #plt.figure()
    #plt.plot(1/li)


    #shield.M_coupling = np.linalg.solve(-shield.inductance, mutual_inductance.T)
    #secondary_C = shield.B_coupling(target_points) @ -shield.M_coupling


    #
    #tmin, tmax = 0.001, 0.001
    #Fs=10000

    time = [0.001, 0.003, 0.005]
    eddy_error = [0.05, 0.01, 0.0025]
    #time_decay = U @ np.exp(-l[None, :]*time[:, None]) @ np.pinv(U)

    time_decay = np.zeros((len(time), shield.inductance.shape[0], shield.inductance.shape[1]))

    induction_spec = []


    Uinv = np.linalg.pinv(U)
    for idx, t in enumerate(time):
         time_decay = U @ np.diag(np.exp(-l*t)) @ Uinv
         eddy_coupling = shield.B_coupling(target_points) @ time_decay @ shield.M_coupling
         induction_spec.append({'coupling':eddy_coupling, 'abs_error':eddy_error[idx], 'rel_error':0, 'target':np.zeros_like(target_field)})





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 4764 vertices by 672 target points... took 1.46 seconds.
    Computing the resistance matrix...




Run QP solver


.. code-block:: default


    import mosek

    coil.s, prob = optimize_streamfunctions(coil,
                                       [target_spec] + induction_spec,
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )

    from bfieldtools.conductor import StreamFunction
    shield.induced_s = StreamFunction(shield.M_coupling @ coil.s, shield)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing the inductance matrix...
    Computing self-inductance matrix using rough quadrature (degree=2). For higher accuracy, set quad_degree to 4 or more.
    Estimating 69923 MiB required for 4764 by 4764 vertices...
    Computing inductance matrix in 180 chunks (8574 MiB memory free), when approx_far=True using more chunks is faster...
    Computing 1/r-potential matrix
    Inductance matrix computation took 127.74 seconds.
    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 16530           
      Cones                  : 1               
      Scalar variables       : 803             
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 16530           
      Cones                  : 1               
      Scalar variables       : 803             
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 401
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 16530             conic                  : 402             
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.25              dense det. time        : 0.00            
    Factor     - ML order time          : 0.00              GP order time          : 0.00            
    Factor     - nonzeros before factor : 8.06e+04          after factor           : 8.06e+04        
    Factor     - dense dim.             : 0                 flops                  : 1.33e+09        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   3.2e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  4.97  
    1   1.9e+01  5.9e-01  1.2e+00  -5.51e-01  6.314966279e+01   6.260479394e+01   5.9e-01  5.08  
    2   1.3e+01  4.0e-01  8.7e-01  -3.18e-01  2.086352147e+02   2.083517407e+02   4.0e-01  5.19  
    3   9.5e+00  2.9e-01  6.2e-01  -1.03e-01  7.110150519e+02   7.108741541e+02   2.9e-01  5.26  
    4   6.9e+00  2.1e-01  5.1e-01  -2.17e-01  8.569354260e+02   8.571030006e+02   2.1e-01  5.36  
    5   2.1e+00  6.5e-02  1.8e-01  -3.63e-01  6.462740473e+03   6.463820814e+03   6.5e-02  5.47  
    6   9.4e-01  2.9e-02  6.2e-02  1.64e-01   1.550188042e+04   1.550256715e+04   2.9e-02  5.55  
    7   4.2e-01  1.3e-02  2.2e-02  3.94e-01   2.266716792e+04   2.266764722e+04   1.3e-02  5.66  
    8   3.6e-01  1.1e-02  1.9e-02  6.82e-01   2.341888345e+04   2.341938164e+04   1.1e-02  5.73  
    9   1.6e-01  5.0e-03  6.7e-03  5.40e-01   2.854405860e+04   2.854438944e+04   5.0e-03  5.86  
    10  4.9e-02  1.5e-03  1.4e-03  6.91e-01   3.274798622e+04   3.274816497e+04   1.5e-03  6.00  
    11  7.8e-03  2.4e-04  1.2e-04  7.38e-01   3.533437395e+04   3.533442513e+04   2.4e-04  6.14  
    12  2.4e-03  7.5e-05  2.1e-05  9.39e-01   3.583090147e+04   3.583091798e+04   7.5e-05  6.23  
    13  9.0e-05  2.8e-06  1.5e-07  9.86e-01   3.605497785e+04   3.605497847e+04   2.8e-06  6.33  
    14  1.3e-06  1.1e-07  9.0e-11  9.99e-01   3.606364888e+04   3.606364888e+04   3.9e-08  6.47  
    15  1.1e-06  8.1e-08  8.1e-11  1.00e+00   3.606366437e+04   3.606366438e+04   3.4e-08  6.70  
    16  1.1e-06  8.1e-08  4.5e-11  1.00e+00   3.606366440e+04   3.606366441e+04   3.4e-08  6.92  
    17  1.1e-06  8.1e-08  4.5e-11  1.00e+00   3.606366440e+04   3.606366441e+04   3.4e-08  7.20  
    18  1.3e-06  4.1e-08  1.4e-11  1.00e+00   3.606371861e+04   3.606371861e+04   1.7e-08  7.41  
    19  1.0e-06  3.5e-08  1.7e-11  1.00e+00   3.606372540e+04   3.606372541e+04   1.5e-08  7.59  
    20  6.7e-07  3.3e-08  2.1e-11  1.00e+00   3.606372838e+04   3.606372838e+04   1.4e-08  7.80  
    21  6.1e-07  3.3e-08  2.3e-11  1.00e+00   3.606372847e+04   3.606372847e+04   1.4e-08  8.02  
    22  4.2e-07  3.1e-08  8.8e-12  1.00e+00   3.606373125e+04   3.606373125e+04   1.3e-08  8.22  
    23  4.2e-07  3.1e-08  1.4e-11  1.00e+00   3.606373125e+04   3.606373126e+04   1.3e-08  8.50  
    24  4.2e-07  3.1e-08  6.2e-12  1.00e+00   3.606373191e+04   3.606373191e+04   1.3e-08  8.70  
    25  5.5e-06  1.5e-08  9.6e-12  1.00e+00   3.606375245e+04   3.606375246e+04   6.4e-09  8.91  
    26  5.4e-06  1.5e-08  1.5e-11  1.00e+00   3.606375249e+04   3.606375250e+04   6.4e-09  9.13  
    27  5.8e-06  7.5e-09  4.0e-12  1.00e+00   3.606376275e+04   3.606376276e+04   3.2e-09  9.30  
    28  5.8e-06  7.5e-09  4.0e-12  1.00e+00   3.606376275e+04   3.606376276e+04   3.2e-09  9.53  
    29  5.8e-06  7.5e-09  4.0e-12  1.00e+00   3.606376275e+04   3.606376276e+04   3.2e-09  9.78  
    30  5.8e-06  7.5e-09  3.1e-12  1.00e+00   3.606376276e+04   3.606376276e+04   3.2e-09  10.05 
    31  2.9e-06  3.8e-09  7.2e-12  1.00e+00   3.606376789e+04   3.606376789e+04   1.6e-09  10.25 
    32  2.9e-06  3.8e-09  7.2e-12  1.00e+00   3.606376789e+04   3.606376789e+04   1.6e-09  10.52 
    33  2.9e-06  3.8e-09  7.2e-12  1.00e+00   3.606376789e+04   3.606376789e+04   1.6e-09  10.75 
    Optimizer terminated. Time: 11.14   


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 3.6063767888e+04    nrm: 7e+04    Viol.  con: 3e-08    var: 0e+00    cones: 0e+00  
      Dual.    obj: 3.6063767889e+04    nrm: 5e+05    Viol.  con: 0e+00    var: 3e-08    cones: 0e+00  




Plot coil windings and target points


.. code-block:: default



    loops, loop_values= scalar_contour(coil.mesh, coil.s.vert, N_contours=6)


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(600, 500))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.005)

    B_target = coil.B_coupling(target_points) @ coil.s

    mlab.quiver3d(*target_points.T, *B_target.T)

    shield.plot_mesh(representation='surface', opacity=0.5, cull_back=True, color=(0.8,0.8,0.8), figure=f)
    shield.plot_mesh(representation='surface', opacity=1, cull_front=True, color=(0.8,0.8,0.8), figure=f)

    f.scene.camera.parallel_projection=1

    f.scene.camera.zoom(1.4)




.. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_004.png
    :class: sphx-glr-single-img





For comparison, let's see how the coils look when we ignore the conducting shield


.. code-block:: default



    coil.unshielded_s, coil.unshielded_prob = optimize_streamfunctions(coil,
                                       [target_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )

    shield.unshielded_induced_s = StreamFunction(shield.M_coupling @ coil.unshielded_s, shield)

    loops, loop_values= scalar_contour(coil.mesh, coil.unshielded_s.vert, N_contours=6)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(600, 500))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.005)

    B_target_unshielded = coil.B_coupling(target_points) @ coil.unshielded_s

    mlab.quiver3d(*target_points.T, *B_target_unshielded.T)

    shield.plot_mesh(representation='surface', opacity=0.5, cull_back=True, color=(0.8,0.8,0.8), figure=f)
    shield.plot_mesh(representation='surface', opacity=1, cull_front=True, color=(0.8,0.8,0.8), figure=f)

    f.scene.camera.parallel_projection=1

    f.scene.camera.zoom(1.4)





.. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_005.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Pre-existing problem not passed, creating...
    Passing parameters to problem...
    Passing problem to solver...


    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 4434            
      Cones                  : 1               
      Scalar variables       : 803             
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer started.
    Problem
      Name                   :                 
      Objective sense        : min             
      Type                   : CONIC (conic optimization problem)
      Constraints            : 4434            
      Cones                  : 1               
      Scalar variables       : 803             
      Matrix variables       : 0               
      Integer variables      : 0               

    Optimizer  - threads                : 8               
    Optimizer  - solved problem         : the dual        
    Optimizer  - Constraints            : 401
    Optimizer  - Cones                  : 1
    Optimizer  - Scalar variables       : 4434              conic                  : 402             
    Optimizer  - Semi-definite variables: 0                 scalarized             : 0               
    Factor     - setup time             : 0.06              dense det. time        : 0.00            
    Factor     - ML order time          : 0.00              GP order time          : 0.00            
    Factor     - nonzeros before factor : 8.06e+04          after factor           : 8.06e+04        
    Factor     - dense dim.             : 0                 flops                  : 3.55e+08        
    ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
    0   3.2e+01  1.0e+00  2.0e+00  0.00e+00   0.000000000e+00   -1.000000000e+00  1.0e+00  1.33  
    1   2.5e+01  7.7e-01  2.5e-01  2.02e+00   3.758490047e+01   3.682062549e+01   7.7e-01  1.36  
    2   1.6e+00  4.9e-02  9.2e-03  1.33e+00   5.286967142e+01   5.284269794e+01   4.9e-02  1.39  
    3   1.0e-01  3.2e-03  1.0e-04  1.07e+00   5.245107002e+01   5.244907373e+01   3.2e-03  1.42  
    4   2.1e-02  6.4e-04  1.1e-05  1.00e+00   5.241313400e+01   5.241274476e+01   6.4e-04  1.45  
    5   1.7e-04  5.3e-06  8.5e-09  1.00e+00   5.241904987e+01   5.241904677e+01   5.3e-06  1.48  
    6   1.7e-06  5.2e-08  8.7e-12  1.00e+00   5.241917851e+01   5.241917848e+01   5.2e-08  1.51  
    7   3.5e-07  6.1e-09  1.5e-13  1.00e+00   5.241918001e+01   5.241917997e+01   1.4e-11  1.55  
    Optimizer terminated. Time: 1.58    


    Interior-point solution summary
      Problem status  : PRIMAL_AND_DUAL_FEASIBLE
      Solution status : OPTIMAL
      Primal.  obj: 5.2419180010e+01    nrm: 1e+02    Viol.  con: 5e-12    var: 0e+00    cones: 0e+00  
      Dual.    obj: 5.2419179969e+01    nrm: 5e+01    Viol.  con: 8e-11    var: 2e-09    cones: 0e+00  




Finally, let's compare the time-courses


.. code-block:: default




    tmin, tmax = 0, 0.025
    Fs=2000

    time = np.linspace(tmin, tmax, int(Fs*(tmax-tmin)+1))

    time_decay = np.zeros((len(time), shield.inductance.shape[0], shield.inductance.shape[1]))

    Uinv = np.linalg.pinv(U)
    for idx, t in enumerate(time):
         time_decay[idx] = U @ np.diag(np.exp(-l*t)) @ Uinv



    B_t = shield.B_coupling(target_points) @ (time_decay @ shield.induced_s).T

    unshieldedB_t = shield.B_coupling(target_points) @ (time_decay @ shield.unshielded_induced_s).T

    import matplotlib.pyplot as plt


    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(8, 4))
    ax.plot(time*1e3, np.mean(np.linalg.norm(B_t, axis=1), axis=0).T, 'k-', label='Minimized', linewidth=1.5)
    ax.set_ylabel('Transient field amplitude')
    ax.semilogy(time*1e3, np.mean(np.linalg.norm(unshieldedB_t, axis=1), axis=0).T, 'k--', label='Ignored', linewidth=1.5 )
    ax.set_xlabel('Time (ms)')


    ax.set_ylim(1e-4, 0.5)
    ax.set_xlim(0, 25)


    plt.grid(which='both', axis='y', alpha=0.1)

    plt.legend()
    fig.tight_layout()

    ax.vlines([1, 5, 10, 20], 1e-4, 0.5, alpha=0.1, linewidth=3, color='r')


.. image:: /auto_examples/coil_design/images/sphx_glr_minimal_eddy_current_cylindrical_coil_design_006.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <matplotlib.collections.LineCollection object at 0x000001D58B0FAC18>




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 8 minutes  20.019 seconds)


.. _sphx_glr_download_auto_examples_coil_design_minimal_eddy_current_cylindrical_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: minimal_eddy_current_cylindrical_coil_design.py <minimal_eddy_current_cylindrical_coil_design.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: minimal_eddy_current_cylindrical_coil_design.ipynb <minimal_eddy_current_cylindrical_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
